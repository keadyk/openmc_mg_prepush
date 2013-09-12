#!/usr/bin/env python
# talliesXML.py

tab = "  " # 2 spaces is one indent-level.

""" This module provides an interface for writing OpenMC input files.  The
classes defined in this module represent the XML input files themselves
(settings, materials, geometry, plots, tallies, and cmfd), and the data
structures owned by these files.  They can be initialized and set up with the 
pertinent data and also handle writing themselves to file. 
"""

valid_filters = ['universe', 'material', 'cell', 'cellborn', 'surface', 'mesh', 
                'energy', 'energyout']

valid_scores = ['flux', 'total', 'scatter', 'nu-scatter', 'scatter-n', 
               'scatter-pn', 'transport', 'diffusion', 'n1n', 'absorption', 
               'fission', 'nu-fission', 'kappa-fission', 'current', 'events']

# These were derived from the ENDF-B/VII data included with the MCNP5 
# RSICC distribution
###TODO:
# - Should probably be replaced with a hash lookup.
# - Also should allow this to be customized based on the user's 
# cross_sections.xml file (and perhaps one set up just for the model at hand).
# - I manually (well, using a built-in function of my text editor) converted
# the nuclide names to lower case.  This should be done automagically.
valid_nuclides = [
  "h-1","h-2","h-3","he-3","he-4","li-6","li-7","be-9","b-10","b-11","c-nat",
  "n-14","n-15","o-16","o-17","f-19","na-22","na-23","mg-24","mg-25","mg-26",
  "al-27","si-28","si-29","si-30","p-31","s-32","s-33","s-34","s-36","cl-35",
  "cl-37","ar-36","ar-38","ar-40","k-39","k-40","k-41","ca-40","ca-42","ca-43",
  "ca-44","ca-46","ca-48","sc-45","ti-46","ti-47","ti-48","ti-49","ti-50","v-nat",
  "cr-50","cr-52","cr-53","cr-54","mn-55","fe-54","fe-56","fe-57","fe-58","co-58",
  "co-58m","co-59","ni-58","ni-59","ni-60","ni-61","ni-62","ni-64","cu-63","cu-65",
  "zn-nat","ga-69","ga-71","ge-70","ge-72","ge-73","ge-74","ge-76","as-74","as-75",
  "se-74","se-76","se-77","se-78","se-79","se-80","se-82","br-79","br-81","kr-78",
  "kr-80","kr-82","kr-83","kr-84","kr-85","kr-86","rb-85","rb-86","rb-87","sr-84",
  "sr-86","sr-87","sr-88","sr-89","sr-90","y-89","y-90","y-91","zr-90","zr-91",
  "zr-92","zr-93","zr-94","zr-95","zr-96","nb-93","nb-94","nb-95","mo-92","mo-94",
  "mo-95","mo-96","mo-97","mo-98","mo-99","mo-100","tc-99","ru-96","ru-98","ru-99",
  "ru-100","ru-101","ru-102","ru-103","ru-104","ru-105","ru-106","rh-103","rh-105",
  "pd-102","pd-104","pd-105","pd-106","pd-107","pd-108","pd-110","ag-107","ag-109",
  "ag-110m","ag-111","cd-106","cd-108","cd-110","cd-111","cd-112","cd-113",
  "cd-114","cd-115m","cd-116","in-113","in-115","sn-112","sn-113","sn-114",
  "sn-115","sn-116","sn-117","sn-118","sn-119","sn-120","sn-122","sn-123",
  "sn-124","sn-125","sn-126","sb-121","sb-123","sb-124","sb-125","sb-126",
  "te-120","te-122","te-123","te-124","te-125","te-126","te-127m","te-128",
  "te-129m","te-130","te-132","i-127","i-129","i-130","i-131","i-135","xe-123",
  "xe-124","xe-126","xe-128","xe-129","xe-130","xe-131","xe-132","xe-133",
  "xe-134","xe-135","xe-136","cs-133","cs-134","cs-135","cs-136","cs-137",
  "ba-130","ba-132","ba-133","ba-134","ba-135","ba-136","ba-137","ba-138",
  "ba-140","la-138","la-139","la-140","ce-136","ce-138","ce-139","ce-140",
  "ce-141","ce-142","ce-143","ce-144","pr-141","pr-142","pr-143","nd-142",
  "nd-143","nd-144","nd-145","nd-146","nd-147","nd-148","nd-150","pm-147",
  "pm-148","pm-148m","pm-149","pm-151","sm-144","sm-147","sm-148","sm-149",
  "sm-150","sm-151","sm-152","sm-153","sm-154","eu-151","eu-152","eu-153",
  "eu-154","eu-155","eu-156","eu-157","gd-152","gd-153","gd-154","gd-155",
  "gd-156","gd-157","gd-158","gd-160","tb-159","tb-160","dy-156","dy-158",
  "dy-160","dy-161","dy-162","dy-163","dy-164","ho-165","ho-166m","er-162",
  "er-164","er-166","er-167","er-168","er-170","lu-175","lu-176","hf-174",
  "hf-176","hf-177","hf-178","hf-179","hf-180","ta-181","ta-182","w-182","w-183",
  "w-184","w-186","re-185","re-187","ir-191","ir-193","au-197","hg-196","hg-198",
  "hg-199","hg-200","hg-201","hg-202","hg-204","pb-204","pb-206","pb-207",
  "pb-208","bi-209","u-232","u-233","u-234","u-235","u-236","u-237","u-238",
  "u-239","u-240","u-241","np-235","np-236","np-237","np-238","np-239","pu-236",
  "pu-237","pu-238","pu-239","pu-240","pu-241","pu-242","pu-243","pu-244",
  "pu-246","ra-223","ra-224","ra-225","ra-226","ac-225","ac-226","ac-227",
  "th-227","th-228","th-229","th-230","th-232","th-233","th-234","pa-231",
  "pa-232","pa-233","am-241","am-242","am-242m","am-243","am-244","am-244m",
  "cm-241","cm-242","cm-243","cm-244","cm-245","cm-246","cm-247","cm-248",
  "cm-249","cm-250","bk-249","bk-250","cf-249","cf-250","cf-251","cf-252",
  "cf-254","es-254","es-255","fm-255","total"]

# These are constants used by the tally class to define the estimator to use.
EST_ANALOG = 0
EST_TRACK = 1
EST_DEFAULT = 2

class tally(object):
  """This class represents a single tally in the tallies.xml file."""
  
  def __init__(self, estimator = EST_DEFAULT, label = None):
    """This function constructs the tally object.  The estimator type and tally 
    label, if desired, should be entered here.
    """
    if ((not isinstance(label, str)) and (label != None)):
      raise TypeError, 'Invalid Label Type!'
    # Set the label
    self.label = label
    
    # Set the estimator. The only valid values are analog, tracklength, or 
    # left as default for the score types.
    if ((estimator == EST_DEFAULT) or (estimator == EST_ANALOG) or 
      (estimator == EST_TRACK)):
      self.estimator = estimator
    else:
      raise ValueError, 'Invalid value for Estimator!'
    
    # Initialize the class values now.
    self.id = -1                  # Will be set during write
    self.filters = []             # Will be appended to by set_filter
    self.bins = []                # Will be appended to by set_filter
    self.nuclides = []            # Will be appended to by set_nuclides 
    self.scores = []              # Will be appended to by set_scores
    
  def set_nuclides(self, nucl):
    """This function sets the tally object's nuclide list. If it is called
    multiple times, the nuclide list will be appended to and not overwritten."""
    
    # Check if a list was provided, since the actions will be different
    # based on a list or a single nuclide.
    if isinstance(nucl, list):
      # If a list is present, check each nuclide for validity before adding.
      for i in xrange(len(nucl)):
        if isinstance(nucl[i], str):
          # Check the nuclide
          nucl[i] = nucl[i].lower()
          if (nucl[i] in valid_nuclides):
            # Add the nuclide to the list for it is valid.
            self.nuclides.append(nucl[i])
          else:
            raise ValueError, 'Invalid Nuclide Entered!'
        else: 
          raise TypeError, 'Invalid Type of Nucl!'
    # Check if a single nuclide was provided as a string
    elif isinstance(nucl, str):
      # Check the nuclide for validity
      nucl = nucl.lower()
      if (nucl in valid_nuclides):
        # Add the nuclide to the list for it is valid.
        self.nuclides.append(nucl)
      else:
        raise ValueError, 'Invalid Nuclide Entered!'
    else:
      raise TypeError, 'Invalid Type of nucl!'
    
  def set_scores(self, scores):
    """This function sets the tally object's score type list. If it is called
    multiple times, the scores list will be appended to and not overwritten."""
    
    # Check if a list was provided, since the actions will be different
    # based on a list or a single score.
    if isinstance(scores, list):
      for i in xrange(len(scores)):
        order = None
        if isinstance(scores[i], str):
          # Check the scores
          scores[i] = scores[i].lower()
          # See if the last one, two, or three digits are numbers.
          # If so, grab the number and replace with an n.
          end = len(scores[i]) - 1
          score_check = scores[i]
          if scores[i][end - 3 : end].isdigit():
            score_check = scores[i][0 : end - 3] + 'n'
          elif scores[i][end - 2 : end].isdigit():
            score_check = scores[i][0 : end - 2] + 'n'
          elif scores[i][end - 1 : end].isdigit():
            score_check = scores[i][0 : end - 1] + 'n'
          if (score_check in valid_scores):
            self.scores.append(scores[i])
          else:
            raise ValueError, 'Invalid Score Entered!'
        else:
          raise TypeError, 'Invalid Type of scores!'
    # Check if a single score was provided as a string
    elif isinstance(scores, str):
      # Check the scores
      scores = scores.lower()
      # See if the last one, two, or three digits are numbers.
      # If so, grab the number and replace with an n.
      end = len(scores) - 1
      score_check = scores
      if scores[end - 3 : end].isdigit():
        score_check = scores[0 : end - 3] + 'n'
      elif scores[end - 2 : end].isdigit():
        score_check = scores[0 : end - 2] + 'n'
      elif scores[end - 1 : end].isdigit():
        score_check = scores[0 : end - 1] + 'n'
      if (score_check in valid_scores):
        self.scores.append(scores)
      else:
        raise ValueError, 'Invalid Score Entered!'
    else:
      raise TypeError, 'Invalid Type of scores!'
    
  def set_filters(self, filters, bins):
    """This function sets the tally object's filter and bin list. If it is 
    called multiple times, the filter and bins list will be appended to and not 
    overwritten."""
    
    if isinstance(bins, list):
      for i in xrange(len(bins)):
        if ((not isinstance(bins[i], float)) and 
          (not isinstance(bins[i], int))):
          raise TypeError, 'Invalid type of bins!'
      if isinstance(filters, list):
        for i in xrange(len(filters)):
          filters[i] = filters[i].lower()
          if (filters[i] in valid_filters):
            self.filters.append(filters[i])
            self.bins.append(bins)
          else:
            pass # ERROR!
      elif isinstance(filters, str):
        filters = filters.lower()
        if (filters in valid_filters):
          self.filters.append(filters)
          self.bins.append(bins)
        else:
          raise ValueError, 'Invalid Filter Entered!'
      else:
        raise TypeError, 'Invalid Type of filters!'
    else:
      raise TypeError, 'Invalid Type of bins!'
      
  def write(self, myFile, id, indent = 1):
    """This function writes this tally object to a given file (myFile)."""
    # Indent is the current level of indentation.
    
    # Set the indentation spaces
    level1 = tab * indent
    level2 = level1 + tab
    level3 = level2 + tab
    level4 = level3 + tab
    
    self.id = id # This is written here since the tallies id
                 # is defined by the order it was written to the file
    myFile.write(level1 + "<tally")
    myFile.write(' id="' + str(self.id).strip() + '">\n')
    
    if self.label != None:
      myFile.write(level2 + '<label>"' + self.label + '"</label>\n')
    
    for f in xrange(len(self.filters)):
      myFile.write(level2 + "<filter>\n" 
        + level3 + "<type>" + self.filters[f].strip() + "</type>\n")
      if len(self.bins[f]) > 0:
        myFile.write(level3 + "<bins>\n")
        for b in xrange(len(self.bins[f])):
          myFile.write(level4 + str(self.bins[f][b]).strip() + "\n")
        myFile.write(level3 + "</bins>\n")
      myFile.write(level2 + "</filter>\n")
      
    myFile.write(level2 + "<scores>\n")
    for s in xrange(len(self.scores)):
      myFile.write(level3 + str(self.scores[s]).strip() + "\n")
    myFile.write(level2 + "</scores>\n")
      
    if self.estimator == EST_ANALOG:
      myFile.write(level2 + "<estimator> analog </estimator>\n")
    elif self.estimator == EST_TRACK:
      myFile.write(level2 + "<estimator> tracklength </estimator>\n")
    
    myFile.write(level1 + "</tally>\n")
    
class assume_sep(object):
  """This class represents the <assume_separate> element of tallies.xml. It 
  has the ability to set initialize itself and write itself to a file in XML
  format."""
  
  def __init__(self, value = False):
    """Construct the object."""
    if isinstance(value, bool):
      self.value = value
    else:
      raise TypeError, 'Invalid Type of value'
    
  def write(self, myFile, indent = 1):
    """Write the <assume_separate> element to the given file in XML format."""
    level1 = tab * indent
    
    myFile.write(level1 + "<assume_separate> " + str(self.value).lower() + 
      " </assume_separate>\n")
 
class mesh(object):
  """This class represents the <mesh> element as used in the tallies.xml file.
  It has the capability to initialize itself and to print itself to a given
  file in the expected XML format.  It is also the base class for mesh types
  needed in settings and cmfd."""
  
  def __init__(self, typeName, lowerLeft, dim, upperRight = None, width = None):
    """This function initializes the mesh, setting all its values according to
    the provided data."""
    # Check for under/over-specified problem
    if (upperRight == None) and (width == None):
      raise ValueError, ('Mesh under-specified, please provide either the ' + 
      '<upper_right> data or the <width> data.')
    if (upperRight != None) and (width != None):
      raise ValueError, ('Mesh under-specified, please provide either the ' + 
      '<upper_right> data or the <width> data.')
    
    self.id = -1
    self.type = typeName
    self.lower_left = None
    self.upper_right = None
    self.width = None
    self.dimension = None
    
    if (len(lowerLeft < 2) or len(lowerleft > 3)):
      self.lower_left = lowerLeft
    else:
      raise ValueError, 'lowerLeft length does not match the required (2 or 3).'
    
    if isinstance(upperRight, list):
      if (len(upperRight < 2) or len(upperRight > 3)):
        self.upper_right = upperRight
      else:
        raise ValueError, 'upperRight length does not match the required (2 or 3).'
    else:
      raise TypeError, 'upperRight is an invalid type.'
      
    if len(dim) == 3:
      self.dimension = dim
    else:
      raise ValueError, 'Dim length does not match the required (3).'
      
    if isinstance(width, list):
      if len(width) == 3:
        self.width = width
      else:
        raise ValueError, 'Width length does not match the required (3).'
  
  def write(self, myFile, id, indent = 1):
    """This function writes the mesh to the given file in the expected XML
    format."""
    level1 = tab * indent
    level2 = level1 + tab
    
    self.id = id
    myFile.write(level1 + '<mesh id="' + str(self.id).strip() + '" type="' + 
      self.type.strip() + '">\n')
    
    myFile.write(level2 + '<dimension>')
    for d in xrange(len(dimension)):
      myFile.write(' ' + str(self.dimension[d]).strip())
    myFile.write('</dimension>\n')
    
    myFile.write(level2 + '<lower_left>')
    for l in xrange(len(lower_left)):
      myFile.write(' ' + str(self.lower_left[l]).strip())
    myFile.write('</lower_left>\n')
    
    if self.upper_right != None:
      myFile.write(level2 + '<upper_right>')
      for u in xrange(len(upper_right)):
        myFile.write(' ' + str(self.upper_right[u]).strip())
      myFile.write('</upper_right>\n')
    
    if self.width != None:
      myFile.write(level2 + '<width>')
      for w in xrange(len(width)):
        myFile.write(' ' + str(self.width[w]).strip())
      myFile.write('</width>\n')
    
    myFile.write('</mesh>\n')

class tallyXML(object):
  """This class represents the complete tallies.xml file and contains all the
  data necessary to write this file."""
  
  def __init__(self, filename, tallies, meshes = None, assume_sep_val = False):
    """Initializes the tallyXML object including setting all the data."""
    
    # Set the filename, do not open until ready to write though.
    self.filename = filename
    # Should error check that these lists that are passed in are all 
    # the right type.
    if isinstance(tallies, tally):
      self.tallies = [tallies]
    elif isinstance(tallies, list):
      self.tallies = tallies
    if isinstance(meshes, mesh):
      self.meshes = [meshes]
    elif isinstance(meshes, list):
      self.meshes = meshes
    elif meshes == None:
      self.meshes = None
      
    self.assume_separate = assume_sep(assume_sep_val)
    
  def write(self):
    """Writes the tallies.xml file to the previously provided filename."""
    myFile = open(self.filename, 'w')
    # Write header
    myFile.write('<?xml version="1.0"?>\n<tallies>\n\n')
    
    # Write Tallies
    for t in xrange(len(self.tallies)):
      self.tallies[t].write(myFile, t + 1)
      myFile.write("\n")
    
    # Write Mesh
    if self.meshes != None:
      for m in xrange(len(self.meshes)):
        self.meshes[i].write(myFile, m)
    
    # Write assume_separate
    self.assume_separate.write(myFile)
    
    # Close File
    myFile.write("\n</tallies>\n")
    myFile.close()
