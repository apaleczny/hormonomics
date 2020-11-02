
#!/usr/bin/python
"""
Copyright (c) 2020 Adam Paleczny

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

==========

Program to parse hormonomics data sets
execute like this: "python hormonomics.py m+h query_file data_set query-mass transform_steps > output_file.csv"
type m or m+h
if dont enter initial compound and number transformation steps will just return a list of matched compounds from the input query file. 

Search for downstream metabolites/synthetic biotransformations where the query file is known mass differences (e.g. lose -OH = 18.010565(M+H dataset)) between the query-mass and the dataset values will be searched, with a tolerance of 0.02 ppm. 
Best match is defined at closest m/z, for transformations, and only the best match is carried forward to the next round. 
Set number of rounds of transformations.
"""
import sys

filter_option = str(sys.argv[1])
query_file = str(sys.argv[2])
dataset_file = str(sys.argv[3])
try:
    compond_subtraction = float(sys.argv[4])
except IndexError:
    compond_subtraction = None
try:
    transform_steps = int(sys.argv[5])
except IndexError:
    transform_steps = 0

tolerance = 0.02 

query_items = list()
query_file = open(query_file, 'r')
for line in query_file.readlines():
    line_split = line.split(',')
    try:
        query_items.append({'name': line_split[0], 
            'm': float(line_split[1]), 'mplush': float(line_split[2])})
    except ValueError:
        pass

if transform_steps < 1:
    mz_index = None
    dataset_file = open(dataset_file, 'r')
    for line in dataset_file.readlines():
        line_split = line.split(',')
        if mz_index == None:
            for i in range(len(line_split)):
                if line_split[i].lower() == 'm.z':
                    mz_index = i
                    print(","+line.strip())
                    break
        else:
            mz_value = float(line_split[mz_index])
            for query_item in query_items:
                if filter_option == 'm':
                    if query_item['m'] - tolerance < mz_value and mz_value < query_item['m'] + tolerance:
                        print(query_item['name'] + "," + line.strip())
                if filter_option == 'm+h':
                    if query_item['mplush'] - tolerance < mz_value and mz_value < query_item['mplush'] + tolerance:
                        print(query_item['name'] + "," + line.strip())
else:
    current_step = 1
    found_compounds = dict()
    mz_index = None
    dataset_file = open(dataset_file, 'r')
    dataset = dataset_file.readlines()
    dataset_file.close()
    for line in dataset:
        line_split = line.split(',')
        if mz_index == None:
            for i in range(len(line_split)):
                if line_split[i].lower() == 'm.z':
                    mz_index = i
                    print(","+line.strip())
                    break
        else:
            mz_value = float(line_split[mz_index])
            if compond_subtraction:
                mz_value = mz_value - compond_subtraction
            for query_item in query_items:
                if filter_option == 'm':
                    search_value = query_item['m']
                elif filter_option == 'm+h':
                    search_value = query_item['mplush']
                if search_value - tolerance < mz_value and mz_value < search_value + tolerance:
                    match_value = abs(search_value - mz_value)
                    try:
                        if match_value < query_item['best_match']:
                            query_item['best_match'] = match_value
                            query_item['best_match_output'] = line.strip()
                    except KeyError:
                        found_compounds[query_item['name']] = {'m': query_item['m'], 'mplush': query_item['mplush'],
                        'best_match': match_value, 'best_match_output': line.strip(), 'cumulative_substitution_value': compond_subtraction*-1}
    print("Transform Step %s" % current_step)
    for key in found_compounds.keys():
        print(key + "," + (found_compounds[key]['best_match_output']))
    current_step += 1

    def find_best_matches(name, substitution_value, cumulative_substitution_value, query_items, dataset, filter_option, tolerance):
        found_compounds = dict()
        mz_index = None
        for line in dataset:
            line_split = line.split(',')
            if mz_index == None:
                for i in range(len(line_split)):
                    if line_split[i].lower() == 'm.z':
                        mz_index = i
                        break
            else:
                mz_value = float(line_split[mz_index])
                if substitution_value:
                    mz_value = mz_value + substitution_value
                mz_value = mz_value + cumulative_substitution_value
                for query_item in query_items:
                    if filter_option == 'm':
                        search_value = query_item['m']
                    elif filter_option == 'm+h':
                        search_value = query_item['mplush']
                    if search_value - tolerance < mz_value and mz_value < search_value + tolerance:
                        match_value = abs(search_value - mz_value)
                        try:
                            if match_value < found_compounds[name+'-'+query_item['name']]['best_match']:
                                found_compounds[name+'-'+query_item['name']]['best_match'] = match_value
                                found_compounds[name+'-'+query_item['name']]['best_match_output'] = line.strip()
                                found_compounds[name+'-'+query_item['name']]['cumulative_substitution_value'] = cumulative_substitution_value + substitution_value + search_value
                        except KeyError:
                            found_compounds[name+'-'+query_item['name']] = {'m': query_item['m'], 'mplush': query_item['mplush'],
                            'best_match': match_value, 'best_match_output': line.strip(),
                            'cumulative_substitution_value': cumulative_substitution_value + substitution_value + search_value}
        return found_compounds
    
    
    new_found_compounts = dict()
    while not (current_step > transform_steps):
        print("Transform Step %s" % current_step)
        for key in found_compounds.keys():
            if filter_option == 'm':
                substitution_value = found_compounds[key]['m']
            elif filter_option == 'm+h':
                substitution_value = found_compounds[key]['mplush']
            cumulative_substitution_value = found_compounds[key]['cumulative_substitution_value']
            new_found_compounds = find_best_matches(key, substitution_value,
             cumulative_substitution_value,
             query_items, dataset, filter_option, tolerance)
        for key in new_found_compounds.keys():
            print(key + "," + (new_found_compounds[key]['best_match_output']))
        current_step += 1
        found_compounds = new_found_compounds
