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

import pandas as pd
import numpy as np    

class dataReaderUS:
    
  def __init__(self, regions, market, renewable_type, start_time, end_time,
                  folder_data):   

    ## different definitions for different markets needed for data reading or 
    ## for class attributes
    if market == 'intraday':
      fcfreq = pd.to_timedelta(15, unit = 'm')  
      keyword_market = 'M15'
      keyword_dt = '15'
      columns_obs = ['type', 'Time', 'power']
      columns_fore = ['type', 'Time']+ ['leadt_0' + x for x in \
        np.char.mod('%d', np.arange(1,10))] + ['leadt_' + x for x in \
        np.char.mod('%d', np.arange(10,21))]
      dropcol = []
    elif market == 'dayahead':
      fcfreq = pd.to_timedelta(24, unit = 'h')
      keyword_market = 'DA'
      keyword_dt = '60'
      columns_obs = ['type', 'Time', 'power', 'aux']
      columns_fore = ['type', 'Time', 'leadt_1', 'aux']
      dropcol = ['aux']
    else:
     raise ValueError('Unrecognized market type: {0}. Valid values are ' + \
       '"intraday" and "dayahead"'.format(market))

    ## different definitions for different renewable types for data reading
    if renewable_type == 'solar':
      renew_type1 = renew_type2 = 'PV'
    elif renewable_type == 'wind':
      renew_type1 = 'WIND'
      renew_type2 = 'Wind'
    else:
      raise ValueError('Unrecognized data type: {0}. Valid values are ' + \
        '"solar" and "wind"'.format(market)) 

    ## read data and create data frames for past and current data
    data = {}
    current_data = {}
    for region in regions:
      ## past data
      print 'load observations:'+region
      filename = folder_data + 'RT_' + renew_type1 + '_' + keyword_dt +\
        'Min_SC0_2022' + '/RT_Region_' + region + '_' + renew_type2 + '_' + \
        keyword_dt + 'Min_SC0_2022_2006.csv'
      aux_data = pd.read_csv(filename, sep=',', header=None)      
      aux_data.columns = columns_obs
      aux_data.index = pd.to_datetime(aux_data.Time)
      aux_data = aux_data['power']
      data[region] = aux_data[start_time:end_time]


      ## current data
      print 'load forecasts:'+region
      filename = folder_data + keyword_market + '_' + renew_type1 + '_' + \
        keyword_dt + 'Min_SC0_2022/'+ keyword_market + '_Region_' + region + \
         '_' + renew_type2 + '_' + keyword_dt + 'Min_SC0_2022_2006.csv'
      aux_data = pd.read_csv(filename, sep=',', header=None)
      aux_data.columns = columns_fore 
      aux_data.index = pd.to_datetime(aux_data.Time)
      aux_data = aux_data.drop(['type', 'Time'], 1)
      aux_data = aux_data.drop(dropcol, axis = 1)
      for ileadt, leadt in enumerate(aux_data.columns, start = 1):
        aux_data_leadt = aux_data[leadt]
        aux_data_leadt.index = aux_data_leadt.index + fcfreq * ileadt
        data[region] = pd.concat([data[region], aux_data_leadt[start_time:end_time]], \
          axis = 1, join = "inner")
      


    ## set object values
    self.data = data  
    self.attributes = {'renewable_type': renewable_type, \
      'start_time': start_time, 'end_time': end_time, 'regions': regions, \
      'market': market, 'fcfreq': fcfreq}
    
    
    
    
     
    
    
    

    
    
    
