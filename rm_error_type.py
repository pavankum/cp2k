#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Vedran Miletic

from __future__ import print_function
from sys import argv
from os import walk
from os.path import join
from re import match, search, sub

__print_debug__ = False

def print_debug(string):
    if __print_debug__:
        print(string)

def find_source_files(directory, exclude_instatiations=False):
    files = []
    template_files = []
    for dirname, dirnames, filenames in walk(directory):
        for filename in filenames:
            if filename[-2:] == '.F' or filename[-4:] == '.f90':
                files.append(join(dirname, filename))
            if filename[-9:] == '.template': # or filename[...:] == '.instantition':
                template_files.append(join(dirname, filename))
    template_names = [template.split('__')[0] for template in template_files]
    filenames_to_remove = []
    if exclude_instatiations:
        for template in template_names:
            for filename in files:
                if template in filename:
                    filenames_to_remove.append(filename)
    for filename in filenames_to_remove:
        files.remove(filename)
    return sorted(files)

def routine_contains_cp_error_type(code_split, i, routine_kind):
    for k in range(1, len(code_split)-i):
        # hack for mincrawl_type used in cp2k/src/swarm/glbopt_mincrawl.F
        if i+k < len(code_split) and ('cp_error_type' in code_split[i+k] or 'mincrawl_type' in code_split[i+k] or 'glbopt_master_type' in code_split[i+k] or 'swarm_master_type' in code_split[i+k] or 'glbopt_worker_type' in code_split[i+k] or 'swarm_worker_type' in code_split[i+k]):
            return True
        elif 'END {0}'.format(routine_kind) in code_split[i+k] or 'END  {0}'.format(routine_kind) in code_split[i+k]:
            break
    return False

def routine_signature_contains_error_variable(code_split, i):
    if 'err' in code_split[i] or  ('&' in code_split[i] and 'err' in code_split[i+1]) or ('&' in code_split[i] and '&' in code_split[i+1] and 'err' in code_split[i+2]) or ('&' in code_split[i] and '&' in code_split[i+1] and '&' in code_split[i+2] and 'err' in code_split[i+3]) or ('&' in code_split[i] and '&' in code_split[i+1] and '&' in code_split[i+2] and '&' in code_split[i+3] and 'err' in code_split[i+4]) or ('&' in code_split[i] and '&' in code_split[i+1] and '&' in code_split[i+2] and '&' in code_split[i+3] and '&' in code_split[i+4] and 'err' in code_split[i+5]) or ('&' in code_split[i] and '&' in code_split[i+1] and '&' in code_split[i+2] and '&' in code_split[i+3] and '&' in code_split[i+4] and '&' in code_split[i+5] and 'err' in code_split[i+6]) or ('&' in code_split[i] and '&' in code_split[i+1] and '&' in code_split[i+2] and '&' in code_split[i+3] and '&' in code_split[i+4] and '&' in code_split[i+5] and '&' in code_split[i+6] and 'err' in code_split[i+7]) or ('&' in code_split[i] and '&' in code_split[i+1] and '&' in code_split[i+2] and '&' in code_split[i+3] and '&' in code_split[i+4] and '&' in code_split[i+5] and '&' in code_split[i+6] and '&' in code_split[i+7] and 'err' in code_split[i+8]) or False:
        return True
    return False

def find_ranges(code_split):
    line_ranges = []
    for i in range(len(code_split)):
        line = code_split[i]
        if '#define' in line or ('!' in line and 'CPPreconditionNoFail' in line):
            if '\\' not in line:
                print_debug("Linija {0} tip {1} sadržaja {2} je pronađena".format(i, '#define', line))
                line_ranges.append((i, i, '#define'))
            else:
                for j in range(i + 1, len(code_split)):
                    end_line = code_split[j]
                    if '\\' not in end_line:
                        print_debug("Linija {0} tip end of {1} sadržaja {2} je pronađena".format(j, '#define', end_line))
                        line_ranges.append((i, j, '#define'))
                        break
                
        number_of_levels = {}
        number_of_levels['SUBROUTINE'] = 0
        number_of_levels['FUNCTION'] = 0
        for method_name in ['SUBROUTINE', 'FUNCTION']:
            if '{0} '.format(method_name) in line and '"{0}'.format(method_name) not in line and '(' in line and routine_contains_cp_error_type(code_split, i, '{0}'.format(method_name)) and 'END {0}' not in line and 'END  {0}' not in line and 'WRITE' not in line:
                # or routine_signature_contains_error_variable(code_split, i)
                # 'WRITE' not in line avoids matching strings containing FUNCTION/SUBROUTINE
                print_debug("Linija {0} tip {1} sadržaja {2} je pronađena".format(i, method_name, line))
                for j in range(i + 1, len(code_split)):
                    end_line = code_split[j]
                    if ('END {0}'.format(method_name) in end_line or 'END  {0}'.format(method_name) in end_line) and number_of_levels[method_name] == 0:
                        print_debug("Linija {0} tip END {1} sadržaja {2} je pronađena".format(j, method_name, end_line))
                        line_ranges.append((i, j, method_name))
                        break

                    for inner_method_name in ['SUBROUTINE', 'FUNCTION']:
                        if '{0} '.format(inner_method_name) in end_line and '"{0}'.format(inner_method_name) not in end_line and '(' in end_line and 'END {0}' not in end_line and 'END  {0}' not in end_line and 'WRITE' not in end_line:
                            print_debug("Linija {0} tip {1} [inner] sadržaja {2} je pronađena".format(j, inner_method_name, end_line))
                            number_of_levels[inner_method_name] += 1
                        elif 'END {0}'.format(inner_method_name) in end_line or 'END  {0}'.format(inner_method_name) in end_line:
                            print_debug("Linija {0} tip END {1} [inner] sadržaja {2} je pronađena".format(j, inner_method_name, end_line))
                            number_of_levels[inner_method_name] -= 1
                            
    print_debug("Raw line ranges are {0}".format(line_ranges))
    return line_ranges

def remove_nested_line_ranges(line_ranges):
    cleaned_line_ranges = line_ranges[:]
    for i in range(len(line_ranges)):
        for j in range(i+1, len(line_ranges)):
            if line_ranges[i][0] < line_ranges[j][0] and line_ranges[i][1] > line_ranges[j][1]:
                print_debug("Line range {1} is a subset of line range {0}".format(line_ranges[i], line_ranges[j]))
                try:
                    cleaned_line_ranges.remove(line_ranges[j])
                except ValueError:
                    print_debug("Line range {0} removed previously".format(line_ranges[j]))
    return cleaned_line_ranges

def clean_line_ranges(line_ranges):
    cleaned_line_ranges = []
    assimilated_line_ranges = []
    for i in range(len(line_ranges)):
        if i < len(line_ranges) - 1 and line_ranges[i] not in assimilated_line_ranges:
            if line_ranges[i][0] < line_ranges[i+1][0] and line_ranges[i+1][1] < line_ranges[i][1]:
                cleaned_line_ranges.append((line_ranges[i][0], line_ranges[i+1][0], line_ranges[i][2]))
                cleaned_line_ranges.append((line_ranges[i+1][0], line_ranges[i+1][1], line_ranges[i][2] + line_ranges[i+1][2]))
                cleaned_line_ranges.append((line_ranges[i+1][1], line_ranges[i][1], line_ranges[i][2]))
                assimilated_line_ranges.append(line_ranges[i+1])
            else:
                cleaned_line_ranges.append(line_ranges[i])
        elif line_ranges[i] not in assimilated_line_ranges:
            cleaned_line_ranges.append(line_ranges[i])
    return cleaned_line_ranges

def remove_error_lines(text_split, line_range):
    unmodified_text_split = text_split[line_range[0]:line_range[1]+1]
    unmodified_text = "\n".join(unmodified_text_split)
    # hacks for CPPreconditionNoFail and mincrawl_type and other Ole's code
    if 'cp_error_type' not in unmodified_text and '#define' not in unmodified_text and ('!' not in unmodified_text and 'CPPreconditionNoFail' not in unmodified_text) and 'mincrawl_type' not in unmodified_text and 'glbopt_master_type' not in unmodified_text and 'swarm_master_type' not in unmodified_text and 'glbopt_worker_type' not in unmodified_text and 'swarm_worker_type' not in unmodified_text:
        print_debug('No cp_error_type and no #define and no !CPPreconditionNoFail nor hacks is in {0}'.format(line_range))
        return unmodified_text_split
    else:
        no_declaration_text_split = []
        error_names = []
        for i in range(len(unmodified_text_split)):
            line = unmodified_text_split[i]
            if ('cp_error_type' in line) or (i > 0 and 'cp_error_type' in unmodified_text_split[i-1] and '&' in unmodified_text_split[i-1]):
                if '&' not in line:
                    error_name = line.split("::")[-1].strip()
                    if ',' in error_name:
                        # multiple error vars declared in same line
                        for single_error_name in error_name.split(','):
                            error_names.append(single_error_name.strip())
                    else:
                        error_names.append(error_name)
                continue
            if 'CALL' in line and '&' not in line and ('cp_error_init' in line or 'cp_error_set' in line or 'cp_error_dealloc_ref' in line or 'cp_error_get' in line):
                continue
            else:
                no_declaration_text_split.append(line)

        # hack for this%error in mincrawl_type
        if 'mincrawl_type' in unmodified_text or 'glbopt_master_type' in unmodified_text:
            error_names.append('this%error')
            #this can be made dynamic with .split("::") if needed
        if 'swarm_master_type' in unmodified_text:
            error_names.append('master%error')
        if 'glbopt_worker_type' in unmodified_text or 'swarm_worker_type' in unmodified_text:
            error_names.append('worker%error')

        no_error_text_split = []
        if not error_names: error_names.append('error')
        print_debug("Using error names: {0}".format(error_names))
        for i in range(len(no_declaration_text_split)):
            line = no_declaration_text_split[i]
            do_append = True
            error_names = sorted(error_names, key=len, reverse=True)
            for error_name in error_names:
                if error_name in line:
                    if line.strip() in ['error={0}'.format(error_name), '{0}=error'.format(error_name), 'error = {0}'.format(error_name), '{0} = error'.format(error_name)]:
                        do_append = False
                        break
                    if 'my_error = force_env' in line or 'worker%error = error' in line:
                        # exact hacks for force_env.F and glbopt_worker.F
                        do_append = False
                        break
                    # nail error arrrays
                    line = sub(r', *{0}\(.*\) *\)'.format(error_name), r')', line)
                    line = sub(r'\( *{0}\(.*\) *,'.format(error_name), r'(', line)
                    line = sub(r'\( *{0}\(.*\) *\)'.format(error_name), r'()', line)
                    line = sub(r', *{0}\(.*\) *,'.format(error_name), r',', line)
                    # nail other kinds of error or error=error, matching on (, ), or ()
                    line = sub(r'\( *({0}|[a-z_]* *= *{0}) *, *'.format(error_name), r'(', line)
                    line = sub(r' *, *({0}|[a-z_]* *= *{0}) *\)'.format(error_name), r')', line)
                    line = sub(r'\( *({0}|[a-z_]* *= *{0}) *\)'.format(error_name), r'()', line)
                    #line = sub(r' * *({0}|[a-z_]* *= *{0}) *\) *, *'.format(error_name), r')', line)
                    line = sub(r' *, *({0}|[a-z_]* *= *{0}) *, *'.format(error_name), r',', line)
                    line = sub(r' ({0}|[a-z_]* *= *{0}) *, *'.format(error_name), r' ', line)
                    if search(', *({0}|[a-z_]* *= *{0})[a-z_=]+ *'.format(error_name), line) is None:
                        line = sub(r', *({0}|[a-z_]* *= *{0}) *'.format(error_name), r'', line)
                    print_debug("Line \"{0}\" has become \"{1}\" after replacements".format(no_declaration_text_split[i], line))
                    if line.strip() in ["&", "ALLOCATE()", "ALLOCATE( )", "ALLOCATE(  )", "DEALLOCATE()", "DEALLOCATE( )", "DEALLOCATE(  )"] or "ALLOCATE({0}".format(error_name) in line:
                        # almost exact hack for ALLOCATE(suberror
                        do_append = False
                        break
                    if match('!\$OMP *&', line.strip()) is not None:
                        # exact hack for '&' equivalent in OpenMP
                        do_append = False
                        break
                    if search('(omp|OMP).*shared\(\)', line) is not None:
                        # clean up OMP error vars properly
                        line = sub(' *shared\( *\) *', '', line)

                    if (search('({0}|[a-z_]* *= *{0}) *\) *'.format(error_name), line) is not None) and '&' in no_declaration_text_split[i-1]:
                        if search('([a-z_]+{0}|[a-z_]* *= *[a-z_]+{0}) *\) *'.format(error_name), line) is not None:
                            # if there is anything other tahn exact matching name, don't do anything
                            continue
                        print_debug("Multiline error situation matched in line {0}".format(line))
                        if no_error_text_split[-1].strip()[0] != '!' or (('!  ' in no_error_text_split[-1].strip()or '!XXX' in no_error_text_split[-1].strip() or '!$' in no_error_text_split[-1].strip()) and '!' in line):
                            # if the line is not commented out
                            if '&' in no_error_text_split[-1]:
                                no_error_text_split[-1] = sub(r', *&', r')', no_error_text_split[-1]).replace(')&', '))').replace(') &', '))').replace('&  )', '))')
                                print_debug("Replacing ampersand in line-1 {0}".format(no_error_text_split[-1]))
                            else:
                                no_error_text_split[-1] += ')'
                        elif len(no_error_text_split) > 1:
                            if no_error_text_split[-2] == '':
                                pass
                            elif (no_error_text_split[-2].strip()[0] != '!' or ('!  ' in no_error_text_split[-2].strip() and '!' in line)) and '&' in no_error_text_split[-2]:
                                no_error_text_split[-2] = sub(r', *&', r')', no_error_text_split[-2]).replace(')&', '))').replace(') &', '))').replace('&  )', '))')
                                print_debug("Replacing ampersand in line-2 {0}".format(no_error_text_split[-2]))
                            elif len(no_error_text_split) > 2:
                                print_debug("Mot replacing ampersand in line-2 {0}".format(no_error_text_split[-2]))
                                if no_error_text_split[-3] == '':
                                    pass
                                elif (no_error_text_split[-3].strip()[0] != '!' or ('!  ' in no_error_text_split[-3].strip() and '!' in line)) and '&' in no_error_text_split[-3]:
                                    no_error_text_split[-3] = sub(r', *&', r')', no_error_text_split[-3]).replace(')&', '))').replace(') &', '))').replace('&  )', '))')
                                    print_debug("Replacing ampersand in line-3 {0}".format(no_error_text_split[-3]))
                                elif len(no_error_text_split) > 6:
                                    # nail one specific case in input_cp2k_thermostats.F
                                    if no_error_text_split[-7] == '':
                                        pass
                                    elif (no_error_text_split[-7].strip()[0] != '!' or ('!  ' in no_error_text_split[-7].strip() and '!' in line)) and '&' in no_error_text_split[-7]:
                                        no_error_text_split[-7] = sub(r', *&', r')', no_error_text_split[-7]).replace(')&', '))').replace(') &', '))').replace('&  )', '))')
                                        print_debug("Replacing ampersand in line-7 {0}".format(no_error_text_split[-7]))
                                    
                                else:
                                    print_debug("Mot replacing ampersand in line-7 {0}".format(no_error_text_split[-7]))
                        line = sub(r'({0}|[a-z_]* *= *{0}) *\) *'.format(error_name), r'', line)
                        if line.strip() in ['', ',&', ', &', ',  &', '!']:
                            do_append = False
                            if ',' in line.strip():
                                no_error_text_split[-1] += ',&'
                            break
                        elif line.strip() in [') EXIT', '* ft', ') THEN', ')THEN', '.LT. &', ') &', ',cp_p_file)', ',cp_p_file)) THEN', '/=0) THEN', '/pint_env%mass(idim)'] or 'RESULT(' in line or 'RESULT (' in line or ',cp_p_file)' in line:
                            # exact hacks
                            no_error_text_split[-1] += '&'
                        if line.strip() in [')']:
                            do_append = False
                            no_error_text_split[-1] += ')'

            if error_name in line:
                print("Possibly unhandled line: {0}".format(line))
            if do_append:
                no_error_text_split.append(line)
        return no_error_text_split

def remove_error_from_text(text, line_ranges = None):
    text_split = text.split('\n')
    if line_ranges is None:
        line_ranges = remove_nested_line_ranges(find_ranges(text_split))

    print_debug("Found line ranges are {0}".format(line_ranges))
    if line_ranges != []:
        text_modified = []
        text_modified += text_split[:line_ranges[0][0]]

        for i in range(len(line_ranges)):
            line_range = line_ranges[i]
            text_modified += remove_error_lines(text_split, line_range)
            if i < len(line_ranges) - 1:
                text_modified += text_split[line_ranges[i][1]+1:line_ranges[i+1][0]]

        text_modified += text_split[line_ranges[-1][1]+1:]
        return "\n".join(text_modified)
    else:
        return text

def remove_error_from_file(a_file, line_ranges):
    descriptor = open(a_file)
    contents = descriptor.read()
    descriptor.close()
    print_debug("Opened file {0}".format(a_file))

    new_contents = remove_error_from_text(contents, line_ranges)
    if new_contents is not None:
        descriptor = open(a_file, "w")
        descriptor.truncate(0)
        descriptor.write(new_contents)
        descriptor.close()
        print_debug("Cleaned up file {0}".format(a_file))

extra_source_files = ['cp2k/src/ewalds_multipole_debug.h',
                      'cp2k/src/semi_empirical_int_args.h',
                      'cp2k/src/semi_empirical_int_debug.h',
                      'cp2k/src/motion/gopt_f77_methods.h',
]
extra_line_ranges = {'cp2k/src/start/cp2k.F': [(38, 322, 'PROGRAM')],
                     'cp2k/src/start/cp2k_shell.F': [(27, 707, 'PROGRAM')],
}

if __name__ == '__main__':
    source_files = []

    try:
        source_files.append(argv[1])
    except IndexError:
        source_files = find_source_files('cp2k/src')
        source_files += extra_source_files
        for a_file in extra_line_ranges.keys():
            remove_error_from_file(a_file, line_ranges=extra_line_ranges[a_file])

    for a_file in source_files:
        remove_error_from_file(a_file, line_ranges=None)
