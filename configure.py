#!/usr/bin/python

import os, sys, re, errno, glob
from collections import defaultdict


# globals
config = {}

modules_to_build = defaultdict(list)

available_modules = []



# define a function that generates the simulation config file
def readConfig(path):  
  execfile(str(path), config)

def fillModuleTree(path):
  for dir in os.listdir(path):
    newpath=path+'/'+dir
    if os.path.isdir(newpath):
      jamfile = glob.glob(path+'/'+dir + '/Jamfile')
      if jamfile:
        fillModuleTree(path+'/'+str(dir))
  available_modules.append(path)

def checkIfModuleIsAvailable(module_name):
  for available_module in available_modules:
    if available_module.endswith('/'+module_name):
      return True
  return False

def getParentPathForModule(module_name):
  for available_module in available_modules:
    if available_module.endswith('/'+module_name):
      return os.path.dirname(available_module)
  
def getJamfileForModule(module_name):
  parent_module_path = str(getParentPathForModule(module_name))
  ending = '/Jamfile'
  if parent_module_path == compwa_root_path:
    ending = '/Jamroot'
  return parent_module_path+ending

def addModuleToBuildDictionary(jamfile_path, module_name):
  submodule_list = modules_to_build[jamfile_path]
  for submodule_name in submodule_list:
    if submodule_name == module_name:
      return
  modules_to_build[jamfile_path].append(module_name)

def recursivelyEnableParentModules(parent_module_path):
  if parent_module_path != compwa_root_path:
    module_name = os.path.basename(parent_module_path)
    parent_dir = os.path.dirname(parent_module_path)
    if checkIfModuleIsAvailable(module_name):
      addModuleToBuildDictionary(getJamfileForModule(module_name), module_name)
      recursivelyEnableParentModules(parent_dir)

def fillModulesToBuildDictionary():
  list_of_modules_to_build = config['requested_modules_to_build']
  for module_name, submodule_list in list_of_modules_to_build.iteritems():
    if submodule_list:
      for submodule in submodule_list:
        relative_path=module_name+"/"+str(submodule)
        if checkIfModuleIsAvailable(relative_path):
          addModuleToBuildDictionary(getJamfileForModule(str(submodule)), str(submodule))
          recursivelyEnableParentModules(os.path.abspath(getParentPathForModule(module_name)+'/'+module_name))
    else:
      if checkIfModuleIsAvailable(module_name):
        addModuleToBuildDictionary(getJamfileForModule(module_name), module_name)    

def readJamfile(url):
  with open(str(url), 'r') as file:
    # read a list of lines into data
    data = file.readlines()
  return data

def writeJamfile(url, data):
  with open(str(url), 'w') as file:
    file.writelines(data)

def enableModule(filedata, module_name):
  searchstring='^.*build-project\s*'+str(module_name)+'\s*;\s*$'
  for index, line in enumerate(filedata):
    result = re.search(searchstring, line)
    if result:
      filedata[index] = re.sub(searchstring, 'build-project '+str(module_name)+' ;\n', line)
    
def disableModule(filedata, module_name):
  searchstring='^\s*build-project\s*'+str(module_name)+'\s*;\s*$'
  for index, line in enumerate(filedata):
    result = re.search(searchstring, line)
    if result:
      filedata[index] = re.sub(searchstring, '#build-project '+str(module_name)+' ;\n', line)
  
def disableAllModules(filedata):
  searchstring='^\s*(build-project\s*\w*\s*;)\s*$'
  for index, line in enumerate(filedata):
    result = re.search(searchstring, line)
    if result:
      filedata[index] = re.sub(searchstring, '#'+result.group(1)+'\n', line)

# main part of code
compwa_root_path=os.path.abspath('.')
jamroot_path=compwa_root_path+'/Jamroot'
if not os.path.isfile(jamroot_path):
  print 'Could not find Jamroot file in current directory! Something is wrong with the project...' 
  exit(1)

# ok at first we should build a tree of modules from the directory structure
fillModuleTree(compwa_root_path)
#print available_modules

# Next step would be to read the config file and make a comparison, but before that:
# The default behaviour is not building a part. Hence everything that is not in the 
# config filem will not be built! It's enought to disable all top level modules
jamroot_data = readJamfile(jamroot_path)
disableAllModules(jamroot_data)
writeJamfile(jamroot_path, jamroot_data)

#read config file 
readConfig('build.cfg')

#check user build request with available modules and add to modules_to_build dict
fillModulesToBuildDictionary()

#loop over dict
for k,v in modules_to_build.iteritems():
  jamfile_data = readJamfile(k)
  disableAllModules(jamfile_data)
  for module_name in v:
    enableModule(jamfile_data, module_name)
  writeJamfile(k, jamfile_data)


