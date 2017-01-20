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

import numpy as np                 
import pandas as pd  
import warnings

from scipy.stats import norm
from scipy.interpolate import interp1d      
from sklearn import linear_model  
import statsmodels.formula.api as smf




## climatology cdf/inv_cdf function 
def clim_cdf(data, loc_NC, max_factor):
  l = len(data)
  probabilities = np.arange(1,l+1)/(float(l)+1)     
  quantiles = np.array(sorted(data))
  quantiles[quantiles < 0] = 0.
  if np.count_nonzero(quantiles)==0:
    #night of solar cases
    quantiles_extended = np.array([loc_NC, 0.0])
    probabilities_extended = np.zeros(len(quantiles_extended))
  else:
    #Extension of quantiles to reach nominal capacity. The climatology 
    #functions are built with observations. This extension is to prevent the cases
    #when forecast are higher than observations and therefore out of range.
    #The value 1.2 is the lowest fit found so far. Could be generalize using directly 
    #the real nominal capacity
    quantiles_extended = np.concatenate([[-1e-5, quantiles.min()/max_factor], \
      quantiles, [quantiles.max()*max_factor, loc_NC*max_factor]])
    probabilities_extended = np.concatenate([[-1e-5,0.],probabilities,[1,1+1e-5]])
  
  clim_cdf = {}
  clim_cdf['cdf'] = interp1d(quantiles_extended, probabilities_extended) 
  clim_cdf['inv_cdf'] = interp1d(probabilities_extended, quantiles_extended)
 
  return clim_cdf


## function to apply climatological cdf to data
def apply_cdf(cdf, data, cdf_keyword):
  tod = np.unique(data.index.time)
  tdata = data.copy()
  for itime in tod:
    tindex = np.array(data.index.time == itime)
    tdata[tindex] = cdf[itime][cdf_keyword](data[tindex])
  return tdata


## cdf from conditional quantile regression
def cqr_cdf(value, quantiles, prob, cdf_keyword): 
  quantiles[quantiles < 0] = 0
  quantiles[quantiles > 1] = 1
  ## extend and sort quantiles
  quantiles[0] = 0
  quantiles[1] = 1
  quantiles = quantiles.sort_index(axis = 1)

  ## extend probabilities
  probabilities = np.concatenate([[0],prob,[1]])
  
  ## compute rank of value in quantiles
  quantiles['value'] = value
  n = quantiles.shape[0]
  rank = quantiles.apply(np.argsort, axis = 1).apply(np.argsort, axis = 1)['value']

  ## get quantiles closest to value with corresponding probabilities
  q = np.array(quantiles)
  x0 = q[range(n), rank-1]
  x1 = q[range(n), rank+1]
  y0 = probabilities[rank-1]
  y1 = probabilities[rank]
  
  ## linear interpolation
  if cdf_keyword == 'cdf':
    interpolation = y0 + (value - x0)*(x1-x0)/(y1- y0)
  elif cdf_keyword == 'inv_cdf':
    interpolation = x0 + (value - y0)*(y1-y0)/(x1- x0)

  return interpolation


## prepare data
def prepare_data(tdata, data, improv_forecast, location):
  leadtimes = data.data[location].columns[1:]
  data2 = {}
  if improv_forecast:
    for ileadt, leadt in enumerate(leadtimes, start = 1):
      ## prepare data
      step = data.attributes['fcfreq'] * ileadt  
      pers = tdata['power'].copy()
      pers.index = pers.index + step
      data2[leadt] = pd.concat( \
        [tdata['power'], tdata[leadt], pers], axis = 1, join = 'inner', \
        keys = ['observations', 'predictions', 'persistence'])
      data2[leadt] = data2[leadt].dropna()
  else:
    for leadt in leadtimes:
      ## prepare data       
      data2[leadt] = pd.concat([tdata['power'], tdata[leadt]], \
        axis = 1, join = 'inner', keys = ['observations', 'predictions'])
      data2[leadt] = data2[leadt].dropna()
  return data2
 

from scenarioGeneration import scenarioGeneration

class modelEstimation:

  def __init__(self, data, gene_cov, improv_forecast): 

    #As observations might not reach nominal capacity of the farms while forecasts 
    #might predict it, necessity to define a factor to multiply to the maximum 
    #of observations in the definition of climatology cdf/inv_cdf
    #This factor will not be needed if NC of farms are included/used
    if data.attributes['renewable_type'] == 'wind': max_factor = 1.05
    elif data.attributes['renewable_type'] == 'solar': max_factor = 1.2

    
    climcdf = {}
    quantreg_mod = {}
    uniform = {}
    improve_mod = {}
    prob = np.concatenate([[0.001],np.arange(0.05,0.951,0.05),[0.999]])

    for location in data.attributes['regions']:
      print 'Fitting model for ' + location
      ## time of day
      tod = np.unique(data.data[location].index.time)

      climcdf[location] = {}
      quantreg_mod[location] = {}
      uniform[location] = {}
      improve_mod[location] = {}
      
     
      ## climatology transform
      print '    Climatology transformation'
      obs = data.data[location]['power']
      loc_NC = max(obs) #close to nominal capacity of farm
      for itime in tod:
        climatology = obs[obs.index.time == itime]
        climcdf[location][itime] = clim_cdf(climatology, loc_NC, max_factor)

      tdata = apply_cdf(cdf = climcdf[location], \
        data = data.data[location], cdf_keyword = 'cdf')


      leadtimes = data.data[location].columns[1:]


      ## prepare data
      data2 = prepare_data(tdata, data, improv_forecast, location)
      ## improve forecast
      if improv_forecast:
        print '    Improve forecast'
        for leadt in leadtimes:
          tindex = data2[leadt].index
          x = data2[leadt][['predictions','persistence']] 
          y = data2[leadt][['observations']] 
          improve_mod[location][leadt] = linear_model.LinearRegression().fit(x,y)
          improved = improve_mod[location][leadt].predict(x)
          data2[leadt]['improved'] = pd.DataFrame(improved, index = tindex)
      else: 
        for leadt in leadtimes:
          data2[leadt]['improved'] = data2[leadt]['predictions']


      ## Quantile regression and data transformation
      print '    Quantile regression'
      for ileadt, leadt in enumerate(leadtimes, start = 1):
        step = data.attributes['fcfreq'] * ileadt
        tindex = data2[leadt].index
        #smf.quantreg generates warning - see documentation for more details
        #warning off just for this section
        warnings.filterwarnings("ignore")
        mod = smf.quantreg('observations ~ improved', data2[leadt])
        quantreg_mod[location][leadt] = {}
        quantiles = pd.DataFrame(index = tindex)
        for q in prob:
          quantreg_mod[location][leadt][q] = mod.fit(q = q)
          quantiles = pd.concat([quantiles, \
            pd.DataFrame({q: quantreg_mod[location][leadt][q].predict()}, \
            index = tindex)], axis = 1)
        warnings.filterwarnings("always")

        ## cdf transform
        obs = data2[leadt]['observations']
        uniform[location][leadt] = cqr_cdf(obs , quantiles, prob, 'cdf')

        ## bring index to issue date
        uniform[location][leadt].index = uniform[location][leadt].index - step


    ## Correlation matrix
    print '    Correlation'
    ## prepare data
    uniform_df = pd.DataFrame({'location': [], 'leadt': [], 'value': []})
    for location in data.attributes['regions']:
      for leadt in uniform[location].keys():
        df_loc_leadT_temp = pd.DataFrame({'location': location, 'leadt': leadt, \
          'value': uniform[location][leadt]})      
        uniform_df = pd.concat([uniform_df, df_loc_leadT_temp])
        del df_loc_leadT_temp
    uniform_pivot = uniform_df.pivot_table(columns=('location','leadt'), \
      values='value')
    norm_df = uniform_df.copy()
    norm_df['value'] = norm.ppf(uniform_df['value'])
    norm_pivot = norm_df.pivot_table(index = norm_df.index, \
      columns =('location','leadt'),values='value')
            
    ## compute correlation matrix       
    correlation_matrix = norm_pivot.corr()
    correlation_matrix[np.isnan(correlation_matrix)] = 0.
    if not np.all(np.diag(correlation_matrix) == 1.):
      print('All diagonal values of correlation matrix are not 1!')
      np.fill_diagonal(correlation_matrix.values, 1.)
            
       
    ## TODO: generalize correlation
    
    self.climcdf = climcdf
    self.improve_mod = improve_mod
    self.quantreg_mod = quantreg_mod
    self.correlation_matrix = correlation_matrix
    self.prob = prob
    self.attributes = {'improv_forecast': improv_forecast, 'gene_cov': gene_cov, \
      'renewable_type': data.attributes['renewable_type'], \
      'fcfreq': data.attributes['fcfreq']}
      
  def scenarios(self, data, nb_scenarios):
    return scenarioGeneration(self, data, nb_scenarios)



