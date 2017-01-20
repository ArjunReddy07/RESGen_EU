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



##INPUTS TO BE MODIFIED
#Path where all .py files are stored
folder_code = 'Code/'
#Path to data
folder_data = 'Data/'
#Output path to store scenarios in csv
folder_output = 'Results/'
#Renewable type to be studied: 'wind' or 'solar'
renewable_type = 'solar' #
#Market type: 'intraday' or 'dayahead'
market = 'intraday'
#regions to be studied
regions = ['SCE','SDGE','PG&E_VLY','PG&E_BAY']
#Starting and ending time of training period ('YYYY-MM-DD HH:MM:SS')
start_time = '2022-01-01 01:00:00'
end_time = '2022-07-01 00:00:00'
#Starting and ending time of testing period - when scenarios will be generated ('YYYY-MM-DD HH:MM:SS')
fore_start_time = '2022-07-19 00:00:00'
fore_end_time = '2022-07-21 00:00:00'
#Use of the improved forecast model (0:no - 1:yes) - only relevant for wind case
improv_forecast = 1 
#Fitting of the covariance matrix through a combination of exponential and Cauchy estimations (0:no - 1:yes)
gene_cov = 0  #TODO: no effect at the moment
#Number of scenarios to be computed
nb_scenarios = 9


##CODE STRUCTURE - DON'T MODIFY IF ONLY USE
import sys
sys.path.insert(0, folder_code)
from dataReaderUS import dataReaderUS
from modelEstimation import modelEstimation
from scenarioGeneration import save_scenarios, save_scenarios_EPRI

data = dataReaderUS(regions, market, renewable_type, start_time, end_time,
  folder_data) 
forec_data = dataReaderUS(regions, market, renewable_type, fore_start_time, 
  fore_end_time, folder_data)             
model = modelEstimation(data, gene_cov, improv_forecast)
scenarios = model.scenarios(forec_data, nb_scenarios)

save_scenarios_EPRI(scenarios, folder_output, IncludeReferenceScenario = True)
