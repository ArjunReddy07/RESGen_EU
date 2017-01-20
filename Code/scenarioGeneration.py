# License: BSD_3_clause
#
# Copyright (c) 2016, Jakob W. Messner, Jan Emil Banning Iversen, 
# Pierre Pinson, Igor Arduin
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.
#
# Neither the name of the Technical University of Denmark (DTU)
# nor the names of its contributors may be used to endorse or
# promote products derived from this software without specific
# prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


from scipy.stats import multivariate_normal, norm
import numpy as np                 
import pandas as pd  
from modelEstimation import apply_cdf, prepare_data
import warnings
import os

# TODO use cqr_cdf from model estimation
from scipy.interpolate import interp1d      

def cqr_cdf(quantiles, prob, cdf_keyword): 
  quantiles[quantiles < 0] = 0
  quantiles[quantiles > 1] = 1
  if np.all(quantiles == 0):
    quantiles_extended = np.concatenate([[0], sorted(quantiles), [0]])
  else:
    quantiles_extended = np.concatenate([[0], sorted(quantiles), [1]])
  probabilities_extended = np.concatenate([[0],prob,[1]])
  if cdf_keyword == 'cdf':
    interpolation = interp1d(quantiles_extended, probabilities_extended)
  elif cdf_keyword == 'inv_cdf':
    interpolation = interp1d(probabilities_extended, quantiles_extended)
  return interpolation

#Function to save scenarios in csv files
def save_scenarios(scenarios,folder_output):
  for idate in scenarios.simulation.keys():
    folder = folder_output+'/' + idate.strftime('%Y%m%d%H%M')
    if not os.path.isdir(folder):
      os.makedirs(folder)
    for iloc in np.unique(scenarios.simulation[idate].keys().droplevel(1)):
      filename = iloc + '_' + idate.strftime('%Y%m%d%H%M')+'.csv'
      scenarios.simulation[idate][iloc].to_csv(folder + '/' + filename)


def save_scenarios_EPRI(scenarios,folder_output, IncludeReferenceScenario = False):
  print 'Save scenarios'
  SCH_TMP = pd.DataFrame()
  SEQ_SCH = pd.DataFrame()
  if scenarios.attributes['renewable_type'] == 'solar':
    renew_type = 'SOLAR'
  elif scenarios.attributes['renewable_type'] == 'wind':
    renew_type = 'WIND'

  for i, idate in enumerate(scenarios.simulation.keys()):
    for iloc in np.unique(scenarios.simulation[idate].keys().droplevel(1)): 
      ScheduleName = [iloc + '_' + x + '_' + idate.strftime('%Y%m%d%H%M') + \
        '_' + renew_type
        for x in np.char.mod('%d', scenarios.simulation[idate][iloc].index[2:])]
      SequenceName = [renew_type + '_' + iloc + '_ST_' + x \
        for x in np.char.mod('%d', scenarios.simulation[idate][iloc].index[2:])]
      if IncludeReferenceScenario:
        ScheduleName.insert(0, iloc + '_' + 'REF' + '_' + idate.strftime('%Y%m%d%H%M'))
        SequenceName.insert(0, renew_type + '_' + iloc + '_ST_REF')
      for ileadt, leadt in enumerate(scenarios.simulation[idate][iloc].columns, start = 1):
        date = idate + ileadt*scenarios.attributes['fcfreq']
        scen = scenarios.simulation[idate][iloc][leadt][2-IncludeReferenceScenario:]
        if IncludeReferenceScenario:
          scen['fc'] = scen[0:].mean()

        SCH_TMP = pd.concat([SCH_TMP, pd.DataFrame({'Schedule': ScheduleName, \
          'Time': date.strftime('%Y.%m.%d %H:%M'), \
          'Value': scen,\
          'Enforce': 1})])
        SEQ_SCH = pd.concat([SEQ_SCH, pd.DataFrame({'Sequence': SequenceName,\
          'Schedule': ScheduleName, 'AvailInt':"", \
          'AvailTime': idate.strftime('%Y.%m.%d %H:%M')})])

  SCH_TMP.to_csv(folder_output + 'SCH_TMP_Stochastic.csv', header = True, sep=',',index=False, 
    columns=['Schedule','Time','Value','Enforce'])
  SEQ_SCH.to_csv(folder_output + 'SEQ_SCH.csv', header = True, sep=',',index=False, 
    columns=['Sequence', 'Schedule', 'AvailInt', 'AvailTime'])



class scenarioGeneration:
  
  def __init__(self, model, data, nb_scenarios):
    print 'Computing scenarios'
    correlation_matrix = model.correlation_matrix
    prob = model.prob
    improv_forecast = model.attributes['improv_forecast']
    n = correlation_matrix.shape[1]
    rv_mvnorm = multivariate_normal(np.zeros(n), correlation_matrix)

    ## multivariate normal simulation 
    ## and prepare data frame to fill final scenarios
    simulation = {}
    simulation2 = {}
    tindex = data.data[data.attributes['regions'][0]].index
    for issue_date in tindex:
      simulation[issue_date] = rv_mvnorm.rvs(nb_scenarios)
      simulation[issue_date] = pd.DataFrame(data = norm.cdf(simulation[issue_date]),
        columns = correlation_matrix.columns)
      simulation2[issue_date] = pd.DataFrame(index = ['obs', 'fc']+ range(nb_scenarios))

    for location in data.attributes['regions']:
      ## climatology transformation 
      tdata = apply_cdf(cdf = model.climcdf[location], \
        data = data.data[location], cdf_keyword = 'cdf')


      leadtimes = data.data[location].columns[1:]

      ## prepare data
      data2 = prepare_data(tdata, data, improv_forecast, location)

      ## improve forecast
      if improv_forecast:
        for leadt in leadtimes:
          tindex3 = data2[leadt].index
          x = data2[leadt][['predictions','persistence']] 
          improved = model.improve_mod[location][leadt].predict(x)
          data2[leadt]['improved'] = pd.DataFrame(improved, index = tindex3)
      else: 
        for leadt in leadtimes:
          data2[leadt]['improved'] = data2[leadt]['predictions']

      ## quantiles
      for ileadt, leadt in enumerate(leadtimes, start = 1):
        tindex2 = data2[leadt].index
        step = data.attributes['fcfreq'] * ileadt

        quantiles = pd.DataFrame(index = tindex2)
        warnings.filterwarnings("ignore")
        for q in prob:
          quantreg_mod = model.quantreg_mod[location][leadt][q]
          quantiletemp = pd.Series(quantreg_mod.predict(exog = data2[leadt]), index = tindex2)
          #If the prediction is 0, we attributes scenarios to be zeros (for solar night)
          quantiletemp[data2[leadt]['predictions'] == 0.0] = 0.0
          quantiles = pd.concat([quantiles, quantiletemp], axis = 1)
        warnings.filterwarnings("always")

        ## cdf transform
        ## TODO: simpler and faster
        for issue_date in tindex:     
          date = issue_date + step   
          if date in tindex:
            obs = [data.data[location]['power'].loc[date]]
            fc = [data.data[location][leadt].loc[date]]
          else:
            obs = [np.nan]
            fc = [np.nan]
          if date in tindex2:
            inv_cdf = cqr_cdf(quantiles.loc[date], prob, 'inv_cdf')
            climcdf_inv = model.climcdf[location][(date).time()]['inv_cdf']
            pred = climcdf_inv(inv_cdf(simulation[issue_date][location][leadt]))
          else:  # fill with NaNs
            pred = pd.Series(np.nan, index = range(nb_scenarios))
            
          ## add obs and fc
          pred2 = pd.Series(np.concatenate([obs, fc, pred]), \
            index = ['obs', 'fc']+ range(nb_scenarios))
          simulation2[issue_date] = pd.concat([simulation2[issue_date], \
            pred2], axis = 1, join = "inner")

    for issue_date in tindex:    
      simulation2[issue_date].columns = correlation_matrix.columns
      

    self.simulation = simulation2
    self.attributes = model.attributes
    

