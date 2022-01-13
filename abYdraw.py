#!/usr/bin/python3
import re
import sys
import os
from PIL import Image
import PIL
import tkinter as tk
from tkinter import filedialog
from tkinter import colorchooser
import tkinter.ttk as ttk


######################################
def Get_input(x):
    input = ""
    with open(x,"r") as f:
        for line in f:
            input += line

    return(input)
######################################
def Get_dictionaries(x):
    """
    takes in IgG in SMILES format and identifies variables for dynamically rendering image
    """

    ###Split chains into dictionaries
    brackets = []
    non_brackets = []
    current_string = ""
    for i in range(len(str(x))):
        if x[i] != "[" and x[i] != "]":
            current_string += str(x[i])
        elif x[i] == "[":
            non_brackets.append(current_string)
            current_string = ""
        elif x[i] == "]":
            brackets.append(current_string)
            current_string = ""
    if x[-1] == "]":
        brackets.append(current_string)
    else:
        non_brackets.append(current_string)


    for i in range(len(non_brackets)):
        non_brackets[i]= re.sub("\s","",non_brackets[i])
    y = ""
    counter = 0
    if len(non_brackets) > len(brackets):
        for i in range(len(brackets)):
            y+=non_brackets[i]
            y+="["
            y+=brackets[i]
            y+="]"
            counter +=1
        y+=non_brackets[counter]
    elif len(brackets) > len(non_brackets):
        for i in range(len(non_brackets)):
            y+=non_brackets[i]
            y+="["
            y+=brackets[i]
            y+="]"
        #y+=non_brackets[counter]
    location_counting = []
    if "|" in y:
        splitx  = y.split("|")
    else:
        splitx = [y]
    if len(splitx) == 4:
        chains  = ['VH.a', 'VL.a', 'VH.b', 'VL.b']
    elif len(splitx) == 3:
        chains  = ['VH.a', 'VL.a','VH.b']
    elif len(splitx) == 2:
        chains  = ['VH.a','VH.b']
    elif len(splitx) == 1:
        chains  = ['VH.a']
    elif len(splitx) ==5:
        chains  = ['VH.a', 'VL.a', 'VH.b', 'VL.b', 'fragment1']
    elif len(splitx) ==6:
        chains  = ['VH.a', 'VL.a', 'VH.b', 'VL.b', 'fragment1', 'fragment2']
    elif len(splitx) ==7:
        chains  = ['VH.a', 'VL.a', 'VH.b', 'VL.b', 'fragment1', 'fragment2', 'fragment3']
    elif len(splitx) ==8:
        chains  = ['VH.a', 'VL.a', 'VH.b', 'VL.b', 'fragment1', 'fragment2', 'fragment3', 'fragment4']

    chain_count = len(splitx)
    VHa     = {}
    VHb     = {}
    VLa     = {}
    VLb     = {}
    disulphide_bridges   = ""
    VHa_VLa_bonds  = ""
    VHb_VLb_bonds  = ""
    CH1a_CL1a_bonds=""
    CH1b_CL1b_bonds=""
    fragment1={}
    fragment2={}
    fragment3={}
    fragment4={}



    for i in range(len(chains)):
        chain = splitx[i].split("-")
        if i < 4:
            if "VH" in chain[0]:
                if len(chain)==1:
                    dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                elif chain[1] != "L":
                    dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                elif chain[1] == "L":
                    if "VH" in chain[2]:
                        dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                    elif "VL" in chain[2]:
                        dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                    else:
                        dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))

            elif "VL" in chain[0]:
                try:
                    if  chain[1] != "L":
                        dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                    elif chain[1] == "L":
                        if len(chain) > 2:
                            if "VH" in chain[2]:
                                dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                            elif "VL" in chain[2]:
                                dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                            else:
                                dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))
                except IndexError:
                    dict = str(re.sub("\.|\+|\_|\*|\!|nano","",str(chains[i])))

            elif "X" in chain[0]:
                 dict = str(re.sub("\.|\+|\_|\*|\!","",str(chains[i])))
            elif "L" in chain[0]:
                dict = str(re.sub("\.|\+|\_|\*|\!","",str(chains[i])))
            else:
                dict = str(re.sub("\.|\+|\_|\*|\!","",str(chains[i])))
        elif i >= 4:
            dict_number = i-3
            dict = str("fragment"+str(dict_number))


        for j in range(len(chain)):
            domain   =  str(chain[j])
            domain = str(re.sub("\[.*\]|\(.*\)|\{.*\}|\[|\'|\]|\.","", str(domain)))
            if domain == "H*":
                domain = re.sub("\*","", str(domain))

            if domain == "L":
                domain = "Linker"
            #print(domain)
            note=""
            if "[" in str(chain[j]) and "]" in str(chain[j]):
                note_list= str(re.findall("\[.*?\]", str(chain[j])))
                note = str(re.sub("\[|\]|\'","", str(note_list)))

            location    = []
            locationstr =  re.findall("\((.*?)\)", str(chain[j]))
            locationstr =  str(re.sub("\[|\'|\]","", str(locationstr)))
            locationstr =  str(re.sub(":",",", str(locationstr)))
            locationstr = list(locationstr.split(","))
            if domain != "Linker":
                if len(locationstr) == 0:
                    error_message = str("ERROR: missing numbering at domain "+domain+domain+str([j])+"\nAll domains must be numbered sequentially from N-terminus to C-terminus")
                    raise_error(lower_canvas, error_message)
                for i in range(len(locationstr)):
                    if i == 0:
                        print(domain, locationstr[i])
                        location_counting.append(int(locationstr[i]))
                    location.append(int(locationstr[i]))

            else:
                location = ['']



            if re.findall("\{.*?\}", str(chain[j])) != []:
                disulphide_bridges = re.findall("\{.*?\}", str(chain[j]))
                disulphide_bridges = re.sub("\{|\'|\}|\[|\]","", str(disulphide_bridges))
                disulphide_bridges = int(disulphide_bridges)
            else:
                disulphide_bridges = 0

            location = [location, disulphide_bridges,note]
            domain = domain+str([j])


            if dict == "VHa" and domain !="":
                VHa[domain] = location
            elif dict == "VHb" and domain !="":
                VHb[domain] = location
            elif dict == "VLa" and domain !="":
                VLa[domain] = location
            elif dict == "VLb" and domain !="":
                VLb[domain] = location
            elif dict == "fragment1" and domain !="":
                fragment1[domain] = location
            elif dict == "fragment2" and domain !="":
                fragment2[domain] = location
            elif dict == "fragment3" and domain !="":
                fragment3[domain] = location
            elif dict == "fragment4" and domain !="":
                fragment4[domain] = location
            else:
                continue
    ###check locations all appear only once and none are missing###



    #print(VHa)
    #print(VLa)
    #print(VHb)
    #print(VLb)
    ###checker###
    VHa_keyslist = list(VHa.keys())
    VLa_keyslist = list(VLa.keys())
    VHb_keyslist = list(VHb.keys())
    VLb_keyslist = list(VLb.keys())
    fragment1_keyslist = list(fragment1.keys())
    fragment2_keyslist = list(fragment2.keys())
    fragment3_keyslist = list(fragment3.keys())
    fragment4_keyslist = list(fragment4.keys())
    VHa_checked = {}
    VHb_checked = {}
    VLa_checked = {}
    VLb_checked = {}
    ##check VH chains interact
    used = []
    chains = [VHa_keyslist,VLa_keyslist,VHb_keyslist,VLb_keyslist,fragment1_keyslist,fragment2_keyslist,fragment3_keyslist,fragment4_keyslist]
    dicts = [VHa,VLa,VHb,VLb,fragment1,fragment2,fragment3,fragment4]
    print(dicts)
    if chain_count > 2:
        checked_heavy_chains = []
        checked_heavy_dicts  = []
        checked_chains_dicts = []
        VHa_VHb_found = False
        for i in range(len(chains)):
            for j in range(len(chains[i])):
                try:
                    current_interactor = dicts[i].get(chains[i][j])[0][1]
                    for a in range(len(chains)):
                        if dicts[a] != dicts[i]:
                            for b in range(len(chains[a])):
                                interactor = dicts[a].get(chains[a][b])[0][0]

                                if current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and VHa_VHb_found == False:
                                    VHa_VHb_found = True
                                    VHa_checked,VHb_checked= dicts[i],dicts[a]
                                    used.append(dicts[i])
                                    used.append(dicts[a])
                                    break
                                elif current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b]))) and VHa_VHb_found == False:
                                    VHa_VHb_found = True
                                    VHa_checked,VHb_checked= dicts[i],dicts[a]
                                    used.append(dicts[i])
                                    used.append(dicts[a])
                                    break
                                elif current_interactor == interactor and ("H[" not in str(chains[i][j]) and "H[" not in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b])))  and VHa_VHb_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])):
                                    VHa_VHb_found = True
                                    VHa_checked,VHb_checked= dicts[i],dicts[a]
                                    used.append(dicts[i])
                                    used.append(dicts[a])
                                    break
                                elif current_interactor == interactor and ("H[" not in str(chains[i]) and "H[" not in str(chains[a])) and (("X[" not in str(chains[i][j]) and "X[" not in str(chains[a][b])) or ("C[" not in str(chains[i][j]) and "C[" not in str(chains[a][b])))  and VHa_VHb_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])) and ("VH" in str(chains[i]) and "VH" in str(chains[a])):
                                    VHa_VHb_found = True
                                    VHa_checked,VHb_checked= dicts[i],dicts[a]
                                    used.append(dicts[i])
                                    used.append(dicts[a])
                                    break
                except IndexError:
                    if ("X" in chains[i][j] or "C[" in chains[i][j]) and len(dicts[i].get(chains[i][j])[0]) == 1 :
                        current_interactor = dicts[i].get(chains[i][j])[0][0]
                        for a in range(len(chains)):
                            if dicts[a] != dicts[i]:
                                for b in range(len(chains[a])):
                                    if ("X" in chains[a][b] or "C" in chains[a][b]):
                                        interactor = dicts[a].get(chains[a][b])[0][0]

                                        if current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and VHa_VHb_found == False:
                                            VHa_VHb_found = True
                                            VHa_checked,VHb_checked= dicts[i],dicts[a]
                                            used.append(dicts[i])
                                            used.append(dicts[a])
                                            break
                                        elif current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b])))  and VHa_VHb_found == False:
                                            VHa_VHb_found = True
                                            VHa_checked,VHb_checked= dicts[i],dicts[a]
                                            used.append(dicts[i])
                                            used.append(dicts[a])
                                            break
                                        elif current_interactor == interactor and ("H[" not in str(chains[i][j]) and "H[" not in str(chains[a][b])) and(("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b]))) and VHa_VHb_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])):
                                            VHa_VHb_found = True
                                            VHa_checked,VHb_checked= dicts[i],dicts[a]
                                            used.append(dicts[i])
                                            used.append(dicts[a])
                                            break
                                        elif current_interactor == interactor and ("H[" not in str(chains[i]) and "H[" not in str(chains[a])) and (("X[" not in str(chains[i][j]) and "X[" not in str(chains[a][b])) or ("C[" not in str(chains[i][j]) and "C[" not in str(chains[a][b])))  and VHa_VHb_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])) and ("VH" in str(chains[i]) and "VH" in str(chains[a])):
                                            VHa_VHb_found = True
                                            VHa_checked,VHb_checked= dicts[i],dicts[a]
                                            used.append(dicts[i])
                                            used.append(dicts[a])
                                            break
                    else:
                        continue

        checked_heavy_dicts.append(VHa_checked)
        checked_heavy_dicts.append(VHb_checked)
        VHa_keyslist = list(VHa_checked.keys())
        VHb_keyslist = list(VHb_checked.keys())
        checked_heavy_chains.append(VHa_keyslist)
        checked_heavy_chains.append(VHb_keyslist)
        VLa_match = False
        VLb_match = False

        for i in range(len(checked_heavy_chains)):
            for j in range(len(checked_heavy_chains[i])):
                try:
                    current_interactor = checked_heavy_dicts[i].get(checked_heavy_chains[i][j])[0][1]
                    for a in range(len(chains)):
                        for b in range(len(chains[a])):
                            if chains[a] not in checked_heavy_chains:
                                interactor = dicts[a].get(chains[a][b])[0][0]
                                print(interactor, current_interactor)
                                if current_interactor == interactor: # and "H[" not in str(chains[i][j]) and "H[" not in str(chains[a]):
                                    if i %2 == 0 and VLa_match == False:
                                        VLa_checked = dicts[a]
                                        VLa_match = True
                                        used.append(dicts[a])
                                    elif i %2 != 0 and VLb_match == False:
                                        VLb_checked = dicts[a]
                                        VLb_match = True
                                        used.append(dicts[a])

                except IndexError:
                    continue
    elif chain_count ==2:
        VHa_checked = VHa
        VHb_checked = VHb
    elif chain_count ==1:
        VHa_checked = VHa
    interacting_fragments = False
    unused = []
    for i in range(len(dicts)):
        found = False
        for j in range(len(used)):
            if dicts[i] == used[j]:
                found = True
                print(dicts[i], used[j])
        if found == False:
            unused.append(dicts[i])
    if chain_count >= 4:
        if VHa_checked == {} or VHb_checked == {} or VLa_checked == {} or VLb_checked == {}:
            error_message = "ERROR: There has been an error in pairing chains in your antibody"
            raise_error(lower_canvas, error_message)
        if fragment1 !={} and fragment2 !={} and fragment3 != {} and fragment4 !={}: #Check for second IgG
        ##Get chains that haven't been used yet
            dicts = [fragment1,fragment2,fragment3,fragment4]
            print("THE UNUSED LOT", unused)
            print("THE USED LOT", used)
            for i in range(len(dicts)):
                for j in range(len(used)):
                    if dicts[i] == used[j]:
                        for k in range(len(unused)):
                            if unused[k] not in dicts:
                                if i==0:
                                    fragment1 = unused[k]
                                elif i==1:
                                    fragment2 = unused[k]
                                elif i==2:
                                    fragment3 = unused[k]
                                elif i==3:
                                    fragment4 = unused[k]

            fragment1_keyslist = list(fragment1.keys())
            fragment2_keyslist = list(fragment2.keys())
            fragment3_keyslist = list(fragment3.keys())
            fragment4_keyslist = list(fragment4.keys())
            fragment1_checked = {}
            fragment2_checked = {}
            fragment3_checked = {}
            fragment4_checked = {}

            ##check VH chains interact
            chains = [fragment1_keyslist,fragment2_keyslist,fragment3_keyslist,fragment4_keyslist]
            dicts = [fragment1,fragment2,fragment3,fragment4]
            if chain_count > 2:
                checked_heavy_chains = []
                checked_heavy_dicts  = []
                checked_chains_dicts = []
                frag1_frag3_found = False
                for i in range(len(chains)):
                    for j in range(len(chains[i])):
                        try:
                            current_interactor = dicts[i].get(chains[i][j])[0][1]
                            for a in range(len(chains)):
                                if dicts[a] != dicts[i]:
                                    for b in range(len(chains[a])):
                                        interactor = dicts[a].get(chains[a][b])[0][0]

                                        if current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and frag1_frag3_found == False:
                                            frag1_frag3_found = True
                                            interacting_fragments = True
                                            fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                            break
                                        elif current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b])))  and frag1_frag3_found == False:
                                            frag1_frag3_found = True
                                            interacting_fragments = True
                                            fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                            break
                                        elif current_interactor == interactor and ("H[" not in str(chains[i][j]) and "H[" not in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b])))  and frag1_frag3_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])):
                                            frag1_frag3_found = True
                                            interacting_fragments = True
                                            fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                            break
                                        elif current_interactor == interactor and ("H[" not in str(chains[i]) and "H[" not in str(chains[a])) and (("X[" not in str(chains[i][j]) and "X[" not in str(chains[a][b])) or ("C[" not in str(chains[i][j]) and "C[" not in str(chains[a][b]))) and frag1_frag3_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])) and ("VH" in str(chains[i]) and "VH" in str(chains[a])):
                                            frag1_frag3_found = True
                                            interacting_fragments = True
                                            fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                            break
                        except IndexError:
                            if ("X" in chains[i][j] or "C[" in chains[i][j]) and len(dicts[i].get(chains[i][j])[0]) == 1 :
                                current_interactor = dicts[i].get(chains[i][j])[0][0]
                                for a in range(len(chains)):
                                    if dicts[a] != dicts[i]:
                                        for b in range(len(chains[a])):
                                            if ("X" in chains[a][b] or "C[" in chains[a][b]):
                                                interactor = dicts[a].get(chains[a][b])[0][0]
                                                if current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and VHa_VHb_found == False:
                                                    frag1_frag3_found = True
                                                    interacting_fragments = True
                                                    fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                                    break
                                                elif current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b])))  and frag1_frag3_found == False:
                                                    frag1_frag3_found = True
                                                    interacting_fragments = True
                                                    fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                                    break
                                                elif current_interactor == interactor and ("H[" not in str(chains[i][j]) and "H[" not in str(chains[a][b])) and (("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) or ("C[" in str(chains[i][j]) and "C[" in str(chains[a][b]))) and frag1_frag3_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])):
                                                    frag1_frag3_found = True
                                                    interacting_fragments = True
                                                    fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                                    break
                                                elif current_interactor == interactor and ("H[" not in str(chains[i]) and "H[" not in str(chains[a])) and (("X[" not in str(chains[i][j]) and "X[" not in str(chains[a][b])) or ("C[" not in str(chains[i][j]) and "C[" not in str(chains[a][b])))  and frag1_frag3_found == False and ("VL" not in str(chains[i]) and "CL" not in str(chains[i]) and "VL" not in str(chains[a]) and  "CL" not in str(chains[a])) and ("VH" in str(chains[i]) and "VH" in str(chains[a])):
                                                    frag1_frag3_found = True
                                                    interacting_fragments = True
                                                    fragment1_checked,fragment3_checked= dicts[i],dicts[a]
                                                    break
                            else:
                                continue
                checked_heavy_dicts.append(fragment1_checked)
                checked_heavy_dicts.append(fragment3_checked)
                fragment1_keyslist = list(fragment1_checked.keys())
                fragment3_keyslist = list(fragment3_checked.keys())
                checked_heavy_chains.append(fragment1_keyslist)
                checked_heavy_chains.append(fragment3_keyslist)
                fragment2_match = False
                fragment4_match = False
                for i in range(len(checked_heavy_chains)):
                    for j in range(len(checked_heavy_chains[i])):
                        try:
                            current_interactor = checked_heavy_dicts[i].get(checked_heavy_chains[i][j])[0][1]
                            for a in range(len(chains)):
                                for b in range(len(chains[a])):
                                    if chains[a] not in checked_heavy_chains:
                                        interactor = dicts[a].get(chains[a][b])[0][0]
                                        if current_interactor == interactor: # and "H[" not in str(chains[i][j]) and "H[" not in str(chains[a]):
                                            if i %2 == 0 and fragment2_match == False:
                                                fragment2_checked = dicts[a]
                                                fragment2_match = True
                                            elif i %2 != 0 and fragment4_match == False:
                                                fragment4_checked = dicts[a]
                                                fragment4_match = True

                        except IndexError:
                            continue

    all_to_check_keys = list(VHa_checked.keys())+list(VLa_checked.keys())+list(VHb_checked.keys())+list(VLb_checked.keys())+list(fragment1.keys())+list(fragment2.keys())+list(fragment3.keys())+list(fragment4.keys())
    for i in range(len(all_to_check_keys)):
        possible_domains = ["VH","VL","CH1","CH2","CH3","CH4","CL","X","H","Linker","L","C"]
        domain_to_print = re.sub("\[.*\]","",all_to_check_keys[i])
        if "V" in str(all_to_check_keys[i]):
            if bool(re.search("[a-h]",all_to_check_keys[i])) == False:
                domain = re.sub("\[.*\]","",all_to_check_keys[i])
                error_message = str("ERROR: Specificity not noted in domain "+domain+"\nAll variable domains must have a specificity attached")
                raise_error(lower_canvas, error_message)
        if "+" in str(all_to_check_keys[i]) and "_" in str(all_to_check_keys[i]):
            error_message= (str("ERROR: Domain cannot be both + and _ "+ domain_to_print))
            raise_error(lower_canvas, error_message)
        if "@" in str(all_to_check_keys[i]) and ">" in str(all_to_check_keys[i]):
            error_message = (str("ERROR: Domain cannot be both @ and > "+domain_to_print))
            raise_error(lower_canvas, error_message)
        if "!" in str(all_to_check_keys[i]) and "CH2" not in str(all_to_check_keys[i]):
            error_message = str("ERROR: ! modifications are only allowed in CH2 domains, not "+ domain_to_print)
            raise_error(lower_canvas, error_message)
        if "Linker[" not in all_to_check_keys[i]:
            domain = re.sub("\.|nano|[a-h]|\@|>|\+|\-|\_|\!|\*","", domain_to_print)
            if domain not in possible_domains:
                print(domain)
                error_message = str("ERROR: Unrecognised domain type "+ str(domain_to_print)+"\nAll domains in expression much be of type VH,VL,CH1,CH2,CH3,CH4,CL,X,H or L")
                raise_error(lower_canvas, error_message)
    if interacting_fragments == False:
        fragment1_checked = fragment1
        fragment2_checked = fragment2
        fragment3_checked = fragment3
        fragment4_checked = fragment4

    if  (fragment1 !={} and fragment2 !={} and fragment3 != {} and fragment4 !={}) and interacting_fragments == True:
        IgG2 = True
    else:
        IgG2 = False
    return(VHa_checked,VLa_checked,VHb_checked,VLb_checked,chain_count,fragment1_checked,fragment2_checked,fragment3_checked,fragment4_checked,IgG2)

######################################
def Check_interactions(chains_list):
#########Set Variables################

    VHa_chain       = chains_list[0]
    VLa_chain       = chains_list[1]
    VHb_chain       = chains_list[2]
    VLb_chain       = chains_list[3]
    print(chains_list[0])
    print(chains_list[1])
    print(chains_list[2])
    print(chains_list[3])
    print(chains_list[5])
    print(chains_list[6])
    print(chains_list[7])
    print(chains_list[8])
    VHa_chain_master       = VHa_chain.copy()
    VLa_chain_master       = VLa_chain.copy()
    VHb_chain_master       = VHb_chain.copy()
    VLb_chain_master       = VLb_chain.copy()
    chain_count     = chains_list[4]
    fragment1       = chains_list[5]
    fragment2       = chains_list[6]
    fragment3       = chains_list[7]
    fragment4       = chains_list[8]
    IgG2            = chains_list[9]
    All_positions_and_chains    ={}
    extra_disulphide_bridges    ={}
    H_disulphide_coordinates    ={}
    completed_disulphidebridges=[]
    Notes           = []
    Notes_positions = []
    noted_Hbonds=False
    H_coordinatey   = []


    def protein_multimer(VHa_chain,VHb_chain,VLa_chain,VLb_chain,fragment1,fragment2,fragment3,fragment4):
    #Get interactors
        domains_in_multimer = {}
        def get_Xs(chain):
            chain_keyslist = list(chain.keys())
            for i in range(len(chain_keyslist)):
                if "X" in chain_keyslist[i]:
                    if len(chain.get(chain_keyslist[i])[0]) > 2:
                        list_of_interactions = []
                        for j in range(len(chain.get(chain_keyslist[i])[0])):
                            list_of_interactions.append(chain.get(chain_keyslist[i])[0][j])
                        domains_in_multimer[chain.get(chain_keyslist[i])[0][0]] = list_of_interactions
        get_Xs(VHa_chain)
        get_Xs(VLa_chain)
        get_Xs(VHb_chain)
        get_Xs(VLb_chain)
        get_Xs(fragment1)
        get_Xs(fragment2)
        get_Xs(fragment3)
        get_Xs(fragment4)

        #check how many multimers are in the antibody
        domains_in_multimer_keyslist = list(domains_in_multimer.keys())

        multimer_combinations = []
        set_multimer_combinations = []
        for i in range(len(domains_in_multimer_keyslist)):
            combination = domains_in_multimer.get(domains_in_multimer_keyslist[i])
            multimer_combinations.append(combination)

        for i in multimer_combinations:
            i.sort()
        set_multimer_combinations = set(tuple(row) for row in multimer_combinations)
        multimer_tuple = list(set_multimer_combinations)

        premade_start_coordinates = {}
        for i in range(len(multimer_tuple)):
            len_multimer = len(multimer_tuple[i])
            for j in range(len(multimer_tuple[i])):
                if j == 0:
                    startx,starty =0,0
                elif j == 1:
                    startx,starty =-30,40
                elif j == 2:
                    startx,starty =60,0
                elif j == 3:
                    startx,starty =-30,40
                elif j == 4:
                    startx,starty = 7,-20
                elif j == 5:
                    startx,starty =-15,0

                premade_start_coordinates[multimer_tuple[i][j]] =[[startx,starty],i+1]

        return(premade_start_coordinates)

    multimers = protein_multimer(VHa_chain,VLa_chain,VHb_chain,VLb_chain,fragment1,fragment2,fragment3,fragment4)
    multimers_keyslist = (multimers.keys())
    finished_multimers = []
    finished_multimers_numbers = []
    finished_multimers_indexes = []
    print(multimers)

    def disulphide_maker(n,bottomx,bottomy,topx,topy,righthanded):
        number_of_bonds = n+1
        height          = topy-bottomy
        if righthanded == True:
            width = bottomx-topx
        elif righthanded == False:
            width = topx-bottomx
        incrementsy      = height/number_of_bonds
        incrementsx      = width/number_of_bonds
        disulphide_points=[]
        currentx = bottomx
        currenty = bottomy
        for i in range(n):
            new_y=currenty+incrementsy
            if righthanded == True:
                new_x=currentx-incrementsx
            elif righthanded == False:
                new_x=currentx+incrementsx
            disulphide_points.append([new_x ,new_y])
            currentx = new_x
            currenty = new_y
        return(disulphide_points)

    def extra_disulphide_maker(n,bottomx,bottomy,topx,topy,righthanded):
        number_of_bonds = n+1
        height          = topy-bottomy
        if righthanded == True:
            width = bottomx-topx
        elif righthanded == False:
            width = topx-bottomx
        incrementsy      = height/number_of_bonds
        incrementsx      = width/number_of_bonds
        disulphide_points=[]
        currentx = bottomx
        currenty = bottomy-10
        for i in range(n):
            if i == 0:
                disulphide_points.append([currentx,currenty])
            elif i > 0:
                new_y=currenty+incrementsy
                if righthanded == True:
                    new_x=currentx+incrementsx
                elif righthanded == False:
                    new_x=currentx-incrementsx
                disulphide_points.append([new_x,new_y])
                currentx = new_x
                currenty = new_y
        return(disulphide_points)

    def assign_to_chain(chain_dict):
    ##Break chain into fragments
        dictionary = chain_dict
        keyslist = list(chain_dict.keys())
        fragments = []
        ADCs      = []
        CCs       = []
        previous_linker = 0

        coordinates_list_heavy_a= []
        coordinates_list_light_a= []
        coordinates_list_heavy_b= []
        coordinates_list_light_b= []
        coordinates_list_heavy_c= []
        coordinates_list_light_c= []
        coordinates_list_heavy_d= []
        coordinates_list_light_d= []
        coordinates_list_heavy_e= []
        coordinates_list_light_e= []
        coordinates_list_heavy_f= []
        coordinates_list_light_f= []
        coordinates_list_heavy_g= []
        coordinates_list_light_g= []
        coordinates_list_heavy_h= []
        coordinates_list_light_h= []
        names_list_heavy_a      = []
        names_list_light_a      = []
        names_list_heavy_b      = []
        names_list_light_b      = []
        names_list_heavy_c      = []
        names_list_light_c      = []
        names_list_heavy_d      = []
        names_list_light_d      = []
        names_list_heavy_e      = []
        names_list_light_e      = []
        names_list_heavy_f      = []
        names_list_light_f      = []
        names_list_heavy_g      = []
        names_list_light_g      = []
        names_list_heavy_h      = []
        names_list_light_h      = []

        for x in range(len(keyslist)):
            print(len(keyslist),x)
            if "Linker" in keyslist[x]:
                fragments.append(keyslist[previous_linker:x])
                previous_linker = x+1
            elif x+1 == len(keyslist):
                fragments.append(keyslist[previous_linker:x+1])
    ##Choose list to append to
        for x in range(len(fragments)):

            current_fragment = []
            for i in range(len(fragments[x])):
                if "H[" not in str(fragments[x][i]) and "X[" not in str(fragments[x][i]):
                    current_fragment.append(fragments[x][i])
                elif "X[" in fragments[x][i]:
                    print("ADC APPENDED")
                    ADCs.append(fragments[x][i])
                elif "C[" in fragments[x][i]:
                    CCs.append(fragments[x][i])
                    print("CC APPENDED")

            fragments[x] = current_fragment
            print(fragments[x])

            def fragment_cleanup(x):
                if "a" in x:
                    x = re.sub("a",".a", x)
                elif "b" in x:
                    x = re.sub("b",".b", x)
                elif "c" in x:
                    x = re.sub("c",".c", x)
                elif "d" in x:
                    x = re.sub("d",".d", x)
                elif "e" in x:
                    x = re.sub("e",".e", x)
                elif "f" in x:
                    x = re.sub("f",".f", x)
                elif "g" in x:
                    x = re.sub("g",".g", x)
                elif "h" in x:
                    x = re.sub("h",".h", x)
                x_cleaned_up = re.sub("\[.*\]","", x)

                return(str(x_cleaned_up))


            if "a" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "d" not in str(fragments[x]) and "e" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_a.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_a.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_a.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_a.append(cleaned_up)
            elif "b" in str(fragments[x]) and "a" not in str(fragments[x]) and "c" not in str(fragments[x]) and "d" not in str(fragments[x]) and "e" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_b.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_b.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_b.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_b.append(cleaned_up)
            elif "c" in str(fragments[x]) and "a" not in str(fragments[x]) and "b" not in str(fragments[x]) and "d" not in str(fragments[x]) and "e" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_c.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_c.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_c.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_c.append(cleaned_up)
            elif "d" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "a" not in str(fragments[x]) and "e" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_d.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_d.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_d.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_d.append(cleaned_up)
            elif "e" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "a" not in str(fragments[x]) and "d" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_e.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_e.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_e.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_e.append(cleaned_up)
            elif "f" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "a" not in str(fragments[x]) and "d" not in str(fragments[x]) and "e" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_f.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_f.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_f.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_f.append(cleaned_up)
            elif "g" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "a" not in str(fragments[x]) and "d" not in str(fragments[x]) and "f" not in str(fragments[x]) and "e" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_g.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_g.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_g.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_g.append(cleaned_up)
            elif "h" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "a" not in str(fragments[x]) and "d" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "e" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_h.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_h.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_h.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_h.append(cleaned_up)
            elif "a" not in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "d" not in str(fragments[x]) and "e" not in str(fragments[x]) and "f" not in str(fragments[x]) and "g" not in str(fragments[x]) and "h" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_a.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_a.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_a.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_a.append(cleaned_up)
            else:
                for i in range(len(fragments[x])):
                    if "a" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_a.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_a.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_a.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_a.append(cleaned_up)
                    elif "b" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_b.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_b.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_b.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_b.append(cleaned_up)
                    elif "c" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_c.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_c.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_c.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_c.append(cleaned_up)
                    elif "d" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_d.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_d.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_d.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_d.append(cleaned_up)
                    elif "e" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_e.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_e.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_e.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_e.append(cleaned_up)
                    elif "f" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_f.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_f.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_f.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_f.append(cleaned_up)
                    elif "g" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_g.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_g.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_g.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_g.append(cleaned_up)
                    elif "h" in str(fragments[x][i]):
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_h.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_h.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_h.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_h.append(cleaned_up)
                    else:
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_a.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_a.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_a.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_a.append(cleaned_up)

        return(coordinates_list_heavy_a,coordinates_list_light_a,coordinates_list_heavy_b,coordinates_list_light_b,coordinates_list_heavy_c,coordinates_list_light_c,coordinates_list_heavy_d,coordinates_list_light_d,names_list_heavy_a,names_list_light_a,names_list_heavy_b,names_list_light_b,names_list_heavy_c,names_list_light_c,names_list_heavy_d,names_list_light_d,coordinates_list_heavy_e,coordinates_list_light_e,coordinates_list_heavy_f,coordinates_list_light_f,coordinates_list_heavy_g,coordinates_list_light_g,coordinates_list_heavy_h,coordinates_list_light_h,names_list_heavy_e,names_list_light_e,names_list_heavy_f,names_list_light_f,names_list_heavy_g,names_list_light_g,names_list_heavy_h,names_list_light_h,CCs)



    def innie_or_outie(chain,VHa_chain,VHb_chain,VLa_chain,VLb_chain,Build_in,Build_out,fragment1,fragment2,fragment3,fragment4, righthanded, IgG2):
        #global IgG2
        innie_or_outie_list = []
        before_H = True
        keyslist = list(chain.keys())
        fragments = fragment1, fragment2, fragment3, fragment4
        fragments_keyslist = [list(fragment1.keys()),list(fragment2.keys()),list(fragment3.keys()),list(fragment4.keys())]
        if "H[" in str(keyslist):
            H_status = True
        else:
            H_status = False


        def in_or_out(n,chain,first_comp,second_comp,Build_in,Build_out,fragment1,fragment2,fragment3,fragment4,Light_chain_check):
            chain_keyslist = list(chain.keys())
            check_VHa_VLa_interactor = False
            if len(chain.get(chain_keyslist[n])[0]) == 1:
                innie_or_outie_list.append("Single_Fv_Chain")
            else:
                interactor = chain.get(chain_keyslist[n])[0][1]
                first_comp_keyslist = list(first_comp.keys())
                for i in range(len(first_comp_keyslist)):
                    current_interactor = first_comp.get(first_comp_keyslist[i])[0][0]
                    if interactor == current_interactor:
                        check_VHa_VLa_interactor = True
                if check_VHa_VLa_interactor == True:
                    innie_or_outie_list.append("outie")
                else:
                    check_VHa_VHb_interactor = False
                    second_comp_keyslist = list(second_comp.keys())
                    for i in range(len(second_comp_keyslist)):
                        current_interactor = second_comp.get(second_comp_keyslist[i])[0][0]
                        if interactor == current_interactor:
                            check_VHa_VHb_interactor = True
                    if check_VHa_VHb_interactor == True:
                        innie_or_outie_list.append("innie")
                    else:
                        check_VH_Fragment_interactor = False
                        for j in range(len(fragments_keyslist)):
                            for f in range(len(fragments_keyslist[j])):
                                if fragments[j] == fragment1:
                                    current_fragment = fragment1
                                elif fragments[j] == fragment2:
                                    current_fragment = fragment2
                                elif fragments[j] == fragment3:
                                    current_fragment = fragment3
                                elif fragments[j] == fragment4:
                                    current_fragment = fragment4
                                print(current_fragment)
                                current_interactor = current_fragment.get(fragments_keyslist[j][f])[0][0]
                                print(current_interactor,interactor)
                                if interactor == current_interactor:
                                    check_VH_Fragment_interactor = True
                        if check_VH_Fragment_interactor == True:
                            innie_or_outie_list.append("outie")
                        elif n > 0:
                            if "Linker[" in keyslist[n-1] or "H[" in keyslist[n-1]:
                                if chain.get(keyslist[n-2])[0][0] == (chain.get(keyslist[n])[0][1]):
                                    if innie_or_outie_list[-2] == "outie":
                                        innie_or_outie_list.append("innie")
                                    elif innie_or_outie_list[-2] == "innie":
                                        innie_or_outie_list.append("outie")
                                elif  chain.get(keyslist[n-2])[0][0] != (chain.get(keyslist[n])[0][1]):
                                    previous_domain_interaction = False
                                    previous_interactor = ""
                                    if "X" not in str(keyslist[n-2]):
                                        interactor = chain.get(keyslist[n-2])[0][1]
                                        for i in range(len(second_comp_keyslist)):
                                            current_interactor = second_comp.get(second_comp_keyslist[i])[0][0]
                                            if interactor == current_interactor:
                                                previous_domain_interaction = True
                                    if (previous_domain_interaction == False and before_H == False) or (previous_domain_interaction == False and H_status == False) or innie_or_outie_list[-2] == "constant":
                                        try:
                                            if chain.get(keyslist[n-2])[0][1] == chain.get(keyslist[n-4])[0][0] and "H[" not in keyslist[n-1]:
                                                if VHb_chain == {}:
                                                    interactor_in_or_out = ""
                                                    interactor = chain.get(keyslist[n])[0][1]
                                                    for i in range(len(keyslist)):
                                                        current_interactor = VHa_chain.get(keyslist[i])[0][0]
                                                        if interactor == current_interactor:
                                                            interactor_in_or_out = innie_or_outie_list[i]
                                                    if interactor_in_or_out == "innie":
                                                        innie_or_outie_list.append("innie")
                                                    elif interactor_in_or_out == "outie":
                                                        innie_or_outie_list.append("innie")
                                                else:
                                                    innie_or_outie_list.append("innie")
                                            else:
                                                innie_or_outie_list.append("outie")
                                        except IndexError:
                                            innie_or_outie_list.append("outie")
                                    elif previous_domain_interaction == True:
                                        if innie_or_outie_list[-2] == "outie":
                                            innie_or_outie_list.append("innie")
                                        elif innie_or_outie_list[-2] == "innie":
                                            innie_or_outie_list.append("outie")
                                    elif innie_or_outie_list[-2] == "outie":
                                        innie_or_outie_list.append("outie")
                                    elif innie_or_outie_list[-2] == "innie":
                                        innie_or_outie_list.append("innie")

        def in_or_out_light(n,chain,first_comp,Build_in,Build_out,Light_chain_check):
            check_VHa_VLa_interactor = False
            chain_keyslist = list(chain.keys())
            interactor = chain.get(chain_keyslist[n])[0][1]
            first_comp_keyslist = list(first_comp.keys())
            for i in range(len(first_comp_keyslist)):
                current_interactor = first_comp.get(first_comp_keyslist[i])[0][0]
                print("OOOARRRGGGG",interactor, current_interactor)
                if interactor == current_interactor:
                    check_VHa_VLa_interactor = True
            if check_VHa_VLa_interactor == True:
                innie_or_outie_list.append("innie")
            elif n > 0:
                if "Linker[" in keyslist[n-1]:
                    if chain.get(keyslist[n-2])[0][0] == (chain.get(keyslist[n])[0][1]):
                        if innie_or_outie_list[-2] == "outie":
                            innie_or_outie_list.append("innie")
                        elif innie_or_outie_list[-2] == "innie":
                            innie_or_outie_list.append("outie")
                    elif  chain.get(keyslist[n-2])[0][0] != (chain.get(keyslist[n])[0][1]):
                        if innie_or_outie_list[-2] == "outie":
                            innie_or_outie_list.append("outie")
                        elif innie_or_outie_list[-2] == "innie":
                            innie_or_outie_list.append("innie")
                        elif innie_or_outie_list[-2] == "constant":
                            innie_or_outie_list.append("outie")


        for n in range(len(chain)):
            print(n, innie_or_outie_list)
            if "V" in keyslist[n]:
                if chain == VHa_chain or chain == VHb_chain:
                    Light_chain_check = False
                    default = "outie"
                elif chain == VLa_chain or chain == VLb_chain:
                    Light_chain_check = True
                    default = "innie"
                elif IgG2 == True and (chain == fragment1 or chain == fragment3):
                    Light_chain_check = False
                    default = "outie"
                else:
                    Light_chain_check = True
                    default = "innie"
                if n == 0:
                    try:
                        if "Linker[" in keyslist[n+1]:
                            if len(chain.get(keyslist[n])[0]) == 1:
                                innie_or_outie_list.append("Single_Fv_Chain")
                            else:
                                try:

                                    if chain.get(keyslist[n])[0][0] == (chain.get(keyslist[n+2])[0][1]) and "Linker[" not in keyslist[n+3]:
                                        if Light_chain_check == False and Build_in == True:
                                            innie_or_outie_list.append("innie")
                                        elif Light_chain_check==False and Build_out== True:
                                            innie_or_outie_list.append("outie")
                                        elif Light_chain_check == True:
                                            innie_or_outie_list.append("innie")
                                    elif chain.get(keyslist[n])[0][0] == (chain.get(keyslist[n+2])[0][1]) and "Linker[" in keyslist[n+3]:
                                        if  Build_out == True:
                                            innie_or_outie_list.append("outie")
                                        elif  Build_in == True:
                                            innie_or_outie_list.append("innie")

                                    elif chain.get(keyslist[n])[0][0] != (chain.get(keyslist[n+2])[0][1]):
                                        if chain_count == 1:
                                            innie_or_outie_list.append("outie")
                                        elif chain_count == 2:
                                            if Light_chain_check == False:
                                                innie_or_outie_list.append("innie")
                                            elif Light_chain_check == True:
                                                innie_or_outie_list.append("innie")
                                        elif chain_count ==4:
                                            innie_or_outie_list.append(default)

                                except IndexError:

                                    innie_or_outie_list.append("innie")

                        elif "Linker[" not in keyslist[n+1]:
                            if len(chain.get(keyslist[n])[0]) == 1:
                                innie_or_outie_list.append("Single_Fv_Chain")
                            elif chain_count == 2 or chain_count == 1:
                                if Light_chain_check == False:
                                    innie_or_outie_list.append("innie")
                            else:
                                innie_or_outie_list.append(default)
                    except IndexError:
                        if len(chain.get(keyslist[n])[0]) == 1:
                            innie_or_outie_list.append("Single_Fv_Chain")
                        else:
                            innie_or_outie_list.append(default)

                elif n > 0:
                    if chain == VHa_chain:
                        in_or_out(n,chain, VLa_chain,VHb_chain,Build_in,Build_out,fragment1,fragment2,fragment3,fragment4,Light_chain_check)
                    elif chain == VLa_chain:
                        in_or_out_light(n,chain,VHa_chain,Build_in,Build_out,Light_chain_check)
                    elif chain == VHb_chain:
                        in_or_out(n,chain, VLb_chain,VHa_chain,Build_in,Build_out,fragment1,fragment2,fragment3,fragment4,Light_chain_check)
                    elif chain == VLb_chain:
                        in_or_out_light(n,chain,VHb_chain,Build_in,Build_out,Light_chain_check)
                    elif IgG2 == True and (chain == fragment1 or chain == fragment3):
                        print("AND THERE YOU HAVE IT FOLKS")
                        innie_or_outie_list.append("outie")
                    else:
                        print("WE DON'T GOT IT FOLKS")
                        innie_or_outie_list.append("innie")


            elif "V" not in keyslist[n] and "H[" not in keyslist[n]:
                innie_or_outie_list.append("constant")
            elif "V" not in keyslist[n] and "H[" in keyslist[n]:
                innie_or_outie_list.append("constant")
                before_H = False
        return(innie_or_outie_list)

    def VH_VL_check(current_chain_pos,current_dictionary,Light_chain_check,Heavy_chain,Light_chain,have,lookingfor):
        interaction = ""
        if Light_chain_check == False:
            keyslist = list(Light_chain.keys())
            candidate_chain = Light_chain
        elif Light_chain_check == True:
            keyslist = list(Heavy_chain.keys())
            candidate_chain = Heavy_chain
        current_chain_loc =  current_dictionary.get(current_chain_pos)[0]
        for x in range(len(keyslist)):
            if lookingfor in keyslist[x]:
                candidate_chain_loc = candidate_chain.get(keyslist[x])[1]
                if current_chain_loc == candidate_chain_loc:
                    interaction = "H_L"
                elif current_chain_loc != candidate_chain_loc:
                    interaction = "H_H"

        return(interaction)

    def find_the_fragment(fragment_start,All_positions_and_chains):
        try:
            start_coordinates = All_positions_and_chains.get(int(fragment_start))
            return(start_coordinates)
        except IndexError:
            pass

    def find_the_fragment2(fragment_start,All_positions_and_chains):
        try:
            start_coordinates = All_positions_and_chains.get(int(fragment_start))
            righthanded = start_coordinates[1]
            if righthanded == True:
                startx=start_coordinates[0][0]+60
                starty=start_coordinates[0][1]
            elif righthanded==False:
                startx=start_coordinates[0][0]-60
                starty=start_coordinates[0][1]
            return([startx,starty])
        except IndexError:
            pass


    def renderchains(dictionary,startx,starty):
        width = lower_canvas.winfo_width()
        height = lower_canvas.winfo_height()
        print(width)
        print(height)
        chain                   = []
        chain_dict              = {}
        coordinates_list_heavy_a= []
        coordinates_list_light_a= []
        coordinates_list_heavy_b= []
        coordinates_list_light_b= []
        coordinates_list_heavy_c= []
        coordinates_list_light_c= []
        coordinates_list_heavy_d= []
        coordinates_list_light_d= []
        names_list_heavy_a      = []
        names_list_light_a      = []
        names_list_heavy_b      = []
        names_list_light_b      = []
        names_list_heavy_c      = []
        names_list_light_c      = []
        names_list_heavy_d      = []
        names_list_light_d      = []
        coordinates_list_heavy_e= []
        coordinates_list_light_e= []
        coordinates_list_heavy_f= []
        coordinates_list_light_f= []
        coordinates_list_heavy_g= []
        coordinates_list_light_g= []
        coordinates_list_heavy_h= []
        coordinates_list_light_h= []
        names_list_heavy_e      = []
        names_list_light_e      = []
        names_list_heavy_f      = []
        names_list_light_f      = []
        names_list_heavy_g      = []
        names_list_light_g      = []
        names_list_heavy_h      = []
        names_list_light_h      = []
        bonds                   = []
        hinges                  = []
        linkers                 = []
        first_interaction = []
        H_disulphide_bridge_count=0
        disulphidebridge1        =[]
        disulphidebridge2        =[]
        disulphidebridge3        =[]
        disulphidebridge4        =[]
        disulphidebridge5        =[]
        Location_Text=[]
        text_coordinates= []
        Domain_Text     = []
        arcs_left       = []
        arcs_right      = []
        arcs_left_slant = []
        arcs_right_slant= []
        Build_up_downlist=[]
        ADCs             =[]
        Chem_con         =[]
        innie_or_outie_list=[]
        interaction_counter = 0
        keyslist = list(dictionary.keys())
        H_count = 0


        if chain_count >= 4 :
            keyslista = list(VHa_chain.keys())
            keyslistb = list(VHb_chain.keys())
            if "H[" not in str(keyslistb) and "H[" not in str(keyslista) and IgG2 == False:
                slant = False
            else:
                slant = True
            if dictionary == VHa_chain or dictionary == VHb_chain:
                keyslist = list(dictionary.keys())
                try:
                    if "Linker[" in keyslist[1] and dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[2])[0][1]):
                        Build_out = True
                        Build_in  = False
                    else:
                        Build_in = True
                        Build_out = False
                except IndexError:
                    Build_in = True
                    Build_out = False
            elif (dictionary == VLb_chain and startx <=600 and starty >=100) or (dictionary == VLa_chain and startx >=150 and starty>=100):
                Build_out = True
                Build_in = False

            else:
                Build_in = True
                Build_out = False

        elif chain_count==3:
            slant = True
            Build_in = True
            Build_out = False
        elif chain_count==2:
            keyslista = list(VHa_chain.keys())
            keyslistb = list(VHb_chain.keys())
            if "H[" in str(dictionary) or "X[" in str(dictionary):
                tangle_found = False
                if chain_count ==2:
                    try:
                        interactor1a = VHa_chain.get(keyslista[0])[0][1]
                        interactor1b = VHa_chain.get(keyslista[0])[0][0]
                        if "Linker[" in keyslista[1]:
                            interactor_2a = VHa_chain.get(keyslista[2])[0][1]
                            interactor_2b = VHa_chain.get(keyslista[2])[0][0]
                        else:
                            interactor_2a = VHa_chain.get(keyslista[1])[0][1]
                            interactor_2b = VHa_chain.get(keyslista[2])[0][0]
                        if interactor_2a < interactor1a and interactor1b != interactor_2a and "H[" in keyslista[-1]:
                            tangle_found = True
                    except IndexError:
                        pass
                if tangle_found == False :
                    if "X" in str(dictionary):
                        if "LEUCINE" in str(dictionary):
                            slant = True
                        else:
                            slant = False
                    elif ("H[" in str(keyslista) and "H[" in str(keyslistb)) or ("H*[" in str(keyslista) and "H*[" in str(keyslistb)) :
                        slant = True
                    else:
                        slant = False
                else:
                    if tangle_found == True:
                        slant = False
                    elif "H[" in str(dictionary) or "H*[" in str(dictionary) :
                        slant = True
                    else:
                        slant = False
            #elif dictionary == VHa_chain and "H[" in str(keyslistb) or dictionary == VHb_chain and  "H[" in str(keyslista):
            #    slant = True
            else:
                slant = False

            keyslista = list(VHa_chain.keys())
            keyslistb = list(VHb_chain.keys())
            try:
                if dictionary == VHb_chain and VHa_chain.get(keyslista[0])[0][0] == (VHb_chain.get(keyslistb[0])[0][1]) and VHa_chain.get(keyslista[0])[0][1] == (VHb_chain.get(keyslistb[0])[0][0]):
                    Build_out = True
                    Build_in = False
                elif dictionary == VHa_chain and VHa_chain.get(keyslista[0])[0][0] == (VHb_chain.get(keyslistb[0])[1]) and VHa_chain.get(keyslista[0])[0][1] == (VHb_chain.get(keyslistb[0])[0]):
                    Build_out = True
                    Build_in = False
                else:
                    print(dictionary)
                    try:
                        if dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[2])[0][1]) and "H[" in keyslist[3]:
                            Build_in = True
                            Build_out = False
                        elif dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[2])[0][1]) and dictionary.get(keyslist[4])[0][0] == (dictionary.get(keyslist[6])[0][1]):
                            Build_out = True
                            Build_in = False
                        elif dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[6])[0][1]) and dictionary.get(keyslist[2])[0][0] == (dictionary.get(keyslist[4])[0][1]):
                            Build_in = False
                            Build_out = True
                        else:
                            Build_in = True
                            Build_out = False
                    except IndexError:
                        Build_in = True
                        Build_out = False
            except IndexError:
                Build_in = False
                Build_out = True


        elif chain_count==1 and  dictionary == VHa_chain:
            slant = False
            Build_out = False
            Build_in = True

        else:
            slant = False
            Build_out = True
            Build_in = False
        before_H = True

        in_out_counter = 0
        righthanded = False
        Light_chain_check = False
        if dictionary == VHb_chain or dictionary == VLb_chain:
            Heavy_chain = VHb_chain
            Light_chain = VLb_chain
            righthanded = True

        elif dictionary == VHa_chain or dictionary == VLa_chain:
            Heavy_chain = VHa_chain
            Light_chain = VLa_chain
            righthanded = False
        if dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2 or dictionary == fragment3 or dictionary == fragment4:
            Light_chain_check = True
        if IgG2 == False:
            if dictionary == fragment1 or dictionary == fragment2 or dictionary == fragment3 or dictionary == fragment4:
                if starty < (height/2):
                    slant = True
                else:
                    slant = False
                if startx < (width/2):
                    righthanded=False
                elif startx>(width/2):
                    righthanded=True
        elif IgG2 == True:
            if dictionary == fragment1 or dictionary == fragment2:
                slant = True
                righthanded = False
                if dictionary == fragment1:
                    Light_chain_check == False
                elif dictionary == fragment2:
                    Light_chain_check == True
            elif dictionary == fragment3 or dictionary == fragment4:
                slant = True
                righthanded = True
                if dictionary == fragment3:
                    Light_chain_check == False
                elif dictionary == fragment4:
                    Light_chain_check == True

        if dictionary == VHa_chain or dictionary == VHb_chain:
            if (len(VHa_chain) == len(VHb_chain)):
                equal_chain_lengths = True
            else:
                equal_chain_lengths = False
        else:
            equal_chain_lengths = True
        innie_or_outie_list = innie_or_outie(dictionary, VHa_chain_master,VHb_chain_master,VLa_chain_master,VLb_chain_master,Build_in,Build_out, fragment1,fragment2,fragment3,fragment4,righthanded, IgG2)
        #print(innie_or_outie_list)

        for i in range(len(dictionary)):
            keyslist = list(dictionary.keys())
            V  = False
            X  = False
            direction = str(innie_or_outie_list[i])
            interaction = ""
            Extra_bond=False
            previous_H =False
            mod= ""
            mod_label=""
            if "V" in keyslist[i]:
                V  = True
            elif "X" in keyslist[i]:
                X  = True
                if "LEUCINE" in str(dictionary.get(keyslist[i])):
                    mod="Leucine"
                elif "humanserumalbumin" in str(dictionary.get(keyslist[i])):
                    mod="HSA"
            elif "C[" in keyslist[i]:
                X = True
                mod = "C"

            Build_up=False
            Build_down=True
            Domain_name = str(re.sub("\@|\>|\<|\[.*\]","",str(keyslist[i])))


            if dictionary.get(keyslist[i])[2] != "":
                note_label = str(dictionary.get(keyslist[i])[2])
                if "TYPE:" not in note_label and "NOTE:" not in note_label and "ANTI:" not in note_label and "LENGTH:" not in note_label and "MOD:" not in note_label:
                    error_message = str("ERROR: Unrecognised comment type given in ["+note_label+"]"+"\nAll comments must start with classifiers TYPE:, MOD:, NOTE:, ANTI: or LENGTH:")
                    raise_error(lower_canvas, error_message)
                Domain_name_to_add = re.sub("\*|\+|\-","",Domain_name)
                Notes.append(Domain_name+" "+note_label)
                if len(Notes_positions) == 0:
                    Notes_positions.append([(width/3),(height-100)])
                elif len(Notes_positions) > 0:
                    XY = (Notes_positions[-1][1])+20
                    Notes_positions.append([(width/3),XY])

            #print(mod_label, mod)

            if dictionary == VHa_chain or dictionary == VHb_chain or dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2 or dictionary == fragment3 or dictionary == fragment4:
                try:
                    location = dictionary.get(keyslist[i])[0][0]
                    disulphide_bridge_count = int(dictionary.get(keyslist[i])[1])
                    location    = dictionary.get(keyslist[i])[0]
                    dictionary[keyslist[i]] = location
                    if "H[" in keyslist[i]:
                        H_disulphide_bridge_count = disulphide_bridge_count

                except:
                    disulphide_bridge_count = 0
                    location = dictionary.get(keyslist[i])[0]
                    dictionary[keyslist[i]] = location
            else:
                disulphide_bridge_count = 0
                location = dictionary.get(keyslist[i])[0]
                dictionary[keyslist[i]] = location
            #print(keyslist[i], location, disulphide_bridge_count)


            if "@" in keyslist[i]:
                mod = "@"
                if Light_chain_check == True or chain_count==2:
                    direction = "outie"
                elif Light_chain_check == False:
                    if before_H == True:
                        direction = "innie"
                    elif before_H==False:
                        direction = "outie"
            elif ">" in keyslist[i]:
                mod = ">"
                if Light_chain_check == True or chain_count==2:
                    direction = "innie"
                elif Light_chain_check == False:
                    if before_H == True:
                        direction = "outie"
                    elif before_H==False:
                        direction = "innie"











            #print(keyslist[i], i, len(dictionary))
            print(keyslist[i])


            if i == 0:
                if "X" in keyslist[i] and dictionary.get(keyslist[i])[0] in multimers_keyslist:
                    current_number = dictionary.get(keyslist[i])[0]
                    if dictionary != VHa_chain and dictionary != VHb_chain and dictionary != VLa_chain and dictionary != VLb_chain and dictionary != fragment1 and dictionary != fragment2:
                        additions = multimers.get(dictionary.get(keyslist[i])[0])[0]
                    else:
                        additions = multimers.get(dictionary.get(keyslist[i])[0])[0]

                    if finished_multimers == []:
                        getcoordinates = domainmaker(All_positions_and_chains,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif additions != [0,0]:
                        multimer_number =  multimers.get(dictionary.get(keyslist[i])[0])[1]
                        previous_index = 0
                        for f in range(len(finished_multimers_numbers)):
                            if finished_multimers_indexes[f] == multimer_number:
                                previous_index = f
                        if righthanded == False:
                            startx = finished_multimers[f][0]-(additions[0]/2)
                        elif righthanded == True:
                            startx = (finished_multimers[f][0]+(additions[0]/2))
                        starty =finished_multimers[f][1]
                        print(startx,starty)
                        getcoordinates = domainmaker(All_positions_and_chains,startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    else:
                        print("EXTRA checkpoint1")
                        getcoordinates = domainmaker(All_positions_and_chains,startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H)

                else:
                    print("checkpoint1")
                    getcoordinates = domainmaker(All_positions_and_chains,startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H)

            elif i > 0:
                previous_domain = keyslist[i-1]
                previous_chain  = chain[i-1]
                if dictionary.get(keyslist[i])[0] != ['']:
                    checker = dictionary.get(keyslist[i])[0]
                else:
                    checker = 0
                if "H[" not in keyslist[i-1] and "Linker[" not in keyslist[i-1]:# and "X" not in keyslist[i-1]:
                    previous_number = (dictionary.get(previous_domain)[0])+1
                elif "X[" in keyslist[i-1] and "Linker[" not in keyslist[i]:
                    print("checkpointX")
                    previous_number = (dictionary.get(keyslist[i])[0])
                    dictionary[keyslist[i-1]] = [previous_number,0]
                    previous_domain = (keyslist[i])
                elif "X" in keyslist[i-1] and chain_count > 2 and "H[" not in str(dictionary) and Light_chain_check == False:
                    slant=False



                if "H[" in keyslist[i-1]:

                    previous_domain = keyslist[i-2]
                    previous_chain  = chain[i-2]
                    previous_number = int(dictionary.get(previous_domain)[0])+2


                if dictionary == VLa_chain or dictionary == VLb_chain:
                    if "CL[" in keyslist[i-1]:
                        slant = False

                #print(keyslist[i], dictionary.get(keyslist[i])[0], previous_number, dictionary.get(keyslist[i])[0] , (dictionary.get(previous_domain)[1]))

                if checker in All_positions_and_chains and ("X" in keyslist[i] or "C[" in keyslist[i]) and len(dictionary.get(keyslist[i])) == 1:
                    startx = All_positions_and_chains.get(checker)[0][0]
                    starty = All_positions_and_chains.get(checker)[0][1]

                    getcoordinates = domainmaker(All_positions_and_chains,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H)

                elif "H[" in keyslist[i]:
                    before_H = False
                    H_count += 1
                    Build_out = True
                    Build_in = False
                    H_coordinatey = bottom_bond

                    if i+1 ==len(dictionary):
                        print("checkpoint2")
                        mod = "H"
                        if dictionary == VHa_chain and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-1])[1]+2):
                            if slant == True and righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant == True and righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant==False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            Extra_bond=True
                            Build_up=True
                            Build_down=False


                        elif dictionary.get(keyslist[i])[1] != (dictionary.get(keyslist[i-1])[1]+2) or dictionary == VHb_chain:
                            print("checkpoint3")


                            if righthanded == True and slant==True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False and slant==True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant==False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)


                elif "X" in keyslist[i] and dictionary.get(keyslist[i])[0] in multimers_keyslist:
                    current_number = dictionary.get(keyslist[i])[0]
                    additions = multimers.get(dictionary.get(keyslist[i])[0])[0]

                    if finished_multimers != []:
                        if righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(finished_multimers[-1][0]-additions[0]),(finished_multimers[-1][1]+additions[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(finished_multimers[-1][0]+additions[0]),(finished_multimers[-1][1]+additions[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif finished_multimers == []:
                        if righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-80),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+80),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                elif "H[" in keyslist[i-1]  and ("X" in keyslist[i] or "C[" in keyslist[i]):
                    print("checkpoint4")
                    if righthanded == False:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif righthanded == True:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)



                elif ("X[" in keyslist[i] or "C[" in keyslist[i]) and "H[" not in keyslist[i-1]:
                    print("checkpoint5")
                    if "Linker[" in keyslist[i-1]:
                        previous_chain = chain[i-2]
                    if "C" in keyslist[i-1] and "C[" not in keyslist[i-1]:
                        if Build_in == True:
                            Build_in = False
                            Build_out = True
                    if chain_count == 2 and  mod !="Leucine":
                        if Build_in == True:
                            Build_in = False
                            Build_out = True
                        if righthanded == False and slant == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and slant == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    elif chain_count !=2 and mod !="Leucine" and "X" not in keyslist[i-1]:
                        if innie_or_outie_list[i-2] == "innie" and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif innie_or_outie_list[i-2] == "outie" and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif innie_or_outie_list[i-2] == "innie" and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif innie_or_outie_list[i-2] == "outie" and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                        elif righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_in = True
                        Build_out = False
                    elif chain_count !=2 and mod !="Leucine" and "X" in keyslist[i-1]:
                        if righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_up = True
                        Build_in = True
                        Build_out = False
                    elif mod =="Leucine":
                        if righthanded == False :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)


                elif "H[" in keyslist[i-1]  and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    print("checkpoint6")
                    previous_H = True
                    slant=False
                    if chain_count <=2:
                        if dictionary.get(keyslist[2]) != [''] and "Linker[" in keyslist[2]:
                            if dictionary.get(keyslist[0])[0] != (dictionary.get(keyslist[2])[1]) :
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                        else:
                            if righthanded == True :
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    elif H_count ==1 :
                        if righthanded == True :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif H_count==2 :
                        if righthanded == True :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    Build_in = False
                    Build_out = True

                elif "Linker[" not in keyslist[i-1] and "X" and len(dictionary.get(keyslist[i-1])) ==1:
                    print("checkpoint7")
                    if "Linker[" in keyslist[i]:
                        pass
                        print("checkpoint7.5")
                    elif slant == True and righthanded == True:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    else:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    Build_in  = False
                    Build_out = True

                elif chain_count == 2 and "H[" not in keyslist[i] and "Linker[" not in keyslist[i-1] and "Linker[" not in keyslist[i] and "X" not in keyslist[i] and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-1])[1]+2):
                    print("checkpoint8")
                    if chain_count == 2:
                        if dictionary == VHa_chain:
                            if slant == True and righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant == True and righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            else:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            Build_up=True
                            Build_down=False
                        elif dictionary == VHb_chain:
                                if slant == True and righthanded == True:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-50,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif slant == True and righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+50,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                else:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    #Build_in  = False
                    #Build_out = True


                elif "H[" not in keyslist[i] and "Linker[" not in keyslist[i-1] and "Linker[" not in keyslist[i] and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    print("checkpoint9")
                    if slant == True and righthanded == True:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    else:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)


                    #Build_in  = False
                    #Build_out = True

##Linker
                elif "Linker[" in keyslist[i-1]:
                    print("checkpoint10")
                    previous_chain = chain[i-2]
                    previous_domain = keyslist[i-2]
##SdFV
                    if len(dictionary.get(keyslist[i])) == 1:
                        print("checkpoint11")
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        else:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_in  = False
                        Build_out = True
                    elif len(dictionary.get(keyslist[i-2])) ==1 and "X" not in keyslist[i-2]:
                        print("checkpoint11")
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        else:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_in  = False
                        Build_out = True

##Self-Interacting chains
                    elif dictionary.get(keyslist[i])[0] == previous_number and str(dictionary.get(keyslist[i])[1]) in Location_Text:
                        print("checkpoint12")
                        to_join_number      = str(dictionary.get(keyslist[i])[1])
                        to_join_coordinates = find_the_fragment(to_join_number,All_positions_and_chains)
                        to_joinx            = to_join_coordinates[0][0]
                        to_joiny            = to_join_coordinates[0][1]
                        to_join_righthanded = to_join_coordinates[1]
                        to_join_direction   = to_join_coordinates[2]
                        if chain_count == 1 and to_joinx <= (width/2) and dictionary.get(keyslist[i])[1] != (dictionary.get(keyslist[i-2])[0]) and len(dictionary) >=8:
                            righthanded = True
                            Build_in = False
                            Build_out = True
                        elif chain_count == 1 and to_joinx >= (width/2)  and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-2])[0]):
                            Build_in = False
                            Build_out = True
                        if to_join_righthanded == False and to_join_direction == 'outie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == False and to_join_direction == 'innie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == False and to_join_direction == 'constant':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == True and to_join_direction == 'outie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == True and to_join_direction == 'innie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == True and to_join_direction == 'constant':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        all_list = list(All_positions_and_chains.keys())
                        get = (All_positions_and_chains.get(all_list[-1]))
                        yprev = get[0][1]


                        if dictionary.get(keyslist[i-2])[0] != dictionary.get(keyslist[i])[1] and to_joiny < yprev:
                            print("checkpoint13")
                            Build_up=True
                            Build_down=False
                        if Build_in == True:
                            print("checkpoint13.1")
                            Build_out =True
                            Build_in = False
                        elif Build_out==True:
                            print("checkpoint13.2")
                            Build_in = True
                            Build_out = False
                        #if change_side == True:
                        #    Build_out =True
                        #    Build_in = False


##ADCs
                    elif ("X" in keyslist[i-2] or "C[" in keyslist[i-2]) :
                        print("checkpoint14")
                        if chain_count ==1:
                            if Build_in == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif righthanded == True:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif Build_out == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif righthanded == True:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif chain_count > 1:
                            if chain_count == 2:
                                slant = False
                            if righthanded == False and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            if chain_count==2:
                                Build_out = True
                                Build_in = False
                            else:
                                Build_in = True
                                Build_out = False
##Build up
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == (dictionary.get(keyslist[0])[1]):
                        print("checkpoint15")
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_up=True
                        Build_down=False
                    elif chain_count == 2 and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and  dictionary.get(keyslist[i])[1]+1 == dictionary.get(keyslist[i-2])[1]:
                        print("checkpoint16")
                        if dictionary == VHa_chain:
                            if slant==True and righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])-45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif  slant==True and righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])+45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant==False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            Build_up=True
                            Build_down=False
                        elif dictionary == VHb_chain:
                            try:
                                if "Linker[" in keyslist[i+1]:
                                    if slant==True and righthanded == True:
                                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])+45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                    elif  slant==True and righthanded == True:
                                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])-45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                    elif slant==False:
                                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                    Build_up=True
                                    Build_down=False
                                else:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            except IndexError:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H)

                        if Build_out == True:
                            Build_in = True
                            Build_out = False
                        elif Build_in == True:
                            Build_out = True
                            Build_in = False

##Build across

                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] == (dictionary.get(previous_domain)[1]):
                        print("checkpoint17")
                        if Build_in == True:
                            if righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif Build_out == True:
                            if righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)

##Build down
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                        print("checkpoint18")
                        if "V" in keyslist[i]:
                            in_out_counter +=1
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        else:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

##H disulphides
            if "H[" in keyslist[i] or "H*[" in keyslist[i]  and (dictionary != VHa_1_test or dictionary != VHb_1_test or dictionary != VLa_1_test or dictionary != VLb_1_test):
                if H_disulphide_bridge_count > 0 and Build_up == False:
                    bottomx=bottom_bond[0]
                    bottomy=bottom_bond[1]
                    if H_count %2 !=0:
                        if slant == True and righthanded == False:
                            topx = bottomx+20
                        elif slant == True and righthanded == True:
                            topx = bottomx-20
                        elif slant==False:
                            topx = bottomx
                        topy = bottomy+40
                    elif H_count %2 == 0:
                        slant = True
                        if slant == True and righthanded == False:
                            topx = bottomx-20
                        elif slant == True and righthanded == True:
                            topx = bottomx+20
                        elif slant==False:
                            topx = bottomx
                        topy = bottomy+40
                    location   = dictionary.get(keyslist[i])[0]
                    interactor = dictionary.get(keyslist[i])[1]
                    H_bond_coordinates = disulphide_maker(H_disulphide_bridge_count,bottomx,bottomy,topx,topy,righthanded)
                    H_disulphide_coordinates[interactor]=H_bond_coordinates
                    H_disulphidebridge_keyslist = list(H_disulphide_coordinates.keys())
                    for j in range(len(H_disulphide_coordinates)):
                        if int(H_disulphidebridge_keyslist[j]) == int(location):
                            for x in range(len(H_disulphide_coordinates.get(interactor))):
                                coordinate1 = H_disulphide_coordinates.get(interactor)[x][0]
                                coordinate2 = H_disulphide_coordinates.get(interactor)[x][1]
                                coordinate3 = H_disulphide_coordinates.get(H_disulphidebridge_keyslist[j])[x][0]
                                coordinate4 = H_disulphide_coordinates.get(H_disulphidebridge_keyslist[j])[x][1]
                                coordinates = [coordinate1,coordinate2,coordinate3,coordinate4]
                                completed_disulphidebridges.append(coordinates)
                    if H_disulphide_bridge_count > 10:
                        if noted_Hbonds==False:
                            Notes.append((H_disulphide_bridge_count, "hinge disulphide bonds"))
                            noted_Hbonds==True
                            if len(Notes_positions) == 0:
                                Notes_positions.append([(width/3),(height-100)])
                            elif len(Notes_positions) > 0:
                                XY = (Notes_positions[-1][1])+20
                                Notes_positions.append([(width/3),XY])






##append coordinates to chains
            if "H[" in keyslist[i] or "Linker[" in keyslist[i] and i+1 == len(keyslist):
                chain.append([])
                bonds.append([])
            #elif "H[" in keyslist[i] and i+1 == len(keyslist):
            #    top_bond = getcoordinates[2]
            #    bonds.append(bottom_bond + top_bond)
            else:
                chain.append(getcoordinates[0])

            if i > 0 and Build_down==True and ("H[" not in keyslist[i] and not "Linker[" in keyslist[i]):
                print("checkpoint19")
                top_bond = getcoordinates[2]
                if "H[" in keyslist[i-1]:
                    hinges.append(bottom_bond + top_bond)
                elif "Linker[" in keyslist[i-1]:
                    linkers.append(bottom_bond + top_bond)
                else:

                    bonds.append(bottom_bond + top_bond)
                if mod=="Leucine":
                    bonds.append(getcoordinates[0])
            elif i > 0 and Build_up==True:
                print("checkpoint20")
                top_bond = getcoordinates[2]
                arc_topx  = top_bond[0]
                arcbottomx= bottom_bond[0]
                if arc_topx != arcbottomx and slant == False:
                    top_bond = getcoordinates[2]
                    arc_topx  = top_bond[0]
                    arc_topy  = top_bond[1]
                    arcbottomx= bottom_bond[0]
                    arcbottomy= bottom_bond[1]
                    if righthanded == True:
                        bonds.append([bottom_bond[0], bottom_bond[1], top_bond[0]-20,top_bond[1]])
                    elif righthanded == False:
                        bonds.append([bottom_bond[0], bottom_bond[1], top_bond[0]+20,top_bond[1]])

                else:
                    if "Linker[" in keyslist[i-1]:
                        Linker= True
                    else:
                        Linker = False
                    if righthanded==False and slant == True:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]
                        arc_topy  = top_bond[1]
                        arcbottomx= bottom_bond[0]+50
                        arcbottomy= bottom_bond[1]
                        arcs_left_slant.append([[arc_topx, arc_topy, arcbottomx,arcbottomy],Linker])
                    elif righthanded==False and slant == False:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]-50
                        arc_topy  = top_bond[1]-20
                        arcbottomx= bottom_bond[0]+50
                        arcbottomy= bottom_bond[1]
                        arcs_left.append([[arc_topx, arc_topy, arcbottomx,arcbottomy],Linker])
                    elif righthanded==True and slant == True:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]
                        arc_topy  = top_bond[1]
                        arcbottomx= bottom_bond[0]-50
                        arcbottomy= bottom_bond[1]
                        arcs_right_slant.append([[arc_topx, arc_topy, arcbottomx,arcbottomy],Linker])
                    elif righthanded==True and slant == False:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]+50
                        arc_topy  = top_bond[1]-20
                        arcbottomx= bottom_bond[0]-50
                        arcbottomy= bottom_bond[1]
                        arcs_right.append([[arc_topx, arc_topy, arcbottomx,arcbottomy],Linker])
                    if Extra_bond==True:
                        extrabondx1=top_bond[0]
                        extrabondy1=top_bond[1]-20
                        extrabondx2=top_bond[0]
                        extrabondy2=top_bond[1]+20
                        extra_bond = [extrabondx1,extrabondy1,extrabondx2,extrabondy2]
                        if "H[" in keyslist[i]:
                            hinges.append(extra_bond)
                        else:
                            bonds.append(extra_bond)
                        if H_disulphide_bridge_count > 0:
                            bottomx=extrabondx1
                            bottomy=extrabondy1
                            topx = extrabondx2
                            topy = extrabondy2
                            location   = dictionary.get(keyslist[i])[0]
                            interactor = dictionary.get(keyslist[i])[1]
                            H_bond_coordinates = disulphide_maker(H_disulphide_bridge_count,bottomx,bottomy,topx,topy,righthanded)
                            H_disulphide_coordinates[interactor]=H_bond_coordinates
                            H_disulphidebridge_keyslist = list(H_disulphide_coordinates.keys())
                            for j in range(len(H_disulphide_coordinates)):
                                if int(H_disulphidebridge_keyslist[j]) == int(location):
                                    for x in range(len(H_disulphide_coordinates.get(interactor))):
                                        coordinate1 = H_disulphide_coordinates.get(interactor)[x][0]
                                        coordinate2 = H_disulphide_coordinates.get(interactor)[x][1]
                                        coordinate3 = H_disulphide_coordinates.get(H_disulphidebridge_keyslist[j])[x][0]
                                        coordinate4 = H_disulphide_coordinates.get(interactor)[x][1]
                                        coordinates = [coordinate1,coordinate2,coordinate3,coordinate4]
                                        completed_disulphidebridges.append(coordinates)









###append domains to dictionary of domain names and coordinates
            chain_dict[keyslist[i]] = getcoordinates[0]

            if i+1 == len(keyslist):
                dictionaries_to_append   = assign_to_chain(chain_dict)
                #print(dictionaries_to_append)
                coordinates_list_heavy_a += dictionaries_to_append[0]
                coordinates_list_light_a += dictionaries_to_append[1]
                coordinates_list_heavy_b += dictionaries_to_append[2]
                coordinates_list_light_b += dictionaries_to_append[3]
                coordinates_list_heavy_c += dictionaries_to_append[4]
                coordinates_list_light_c += dictionaries_to_append[5]
                coordinates_list_heavy_d += dictionaries_to_append[6]
                coordinates_list_light_d += dictionaries_to_append[7]
                names_list_heavy_a       += dictionaries_to_append[8]
                names_list_light_a       += dictionaries_to_append[9]
                names_list_heavy_b       += dictionaries_to_append[10]
                names_list_light_b       += dictionaries_to_append[11]
                names_list_heavy_c       += dictionaries_to_append[12]
                names_list_light_c       += dictionaries_to_append[13]
                names_list_heavy_d       += dictionaries_to_append[14]
                names_list_light_d       += dictionaries_to_append[15]
                coordinates_list_heavy_e += dictionaries_to_append[16]
                coordinates_list_light_e += dictionaries_to_append[17]
                coordinates_list_heavy_f += dictionaries_to_append[18]
                coordinates_list_light_f += dictionaries_to_append[19]
                coordinates_list_heavy_g += dictionaries_to_append[20]
                coordinates_list_light_g += dictionaries_to_append[21]
                coordinates_list_heavy_h += dictionaries_to_append[22]
                coordinates_list_light_h += dictionaries_to_append[23]
                names_list_heavy_e       += dictionaries_to_append[24]
                names_list_light_e       += dictionaries_to_append[25]
                names_list_heavy_f       += dictionaries_to_append[26]
                names_list_light_f       += dictionaries_to_append[27]
                names_list_heavy_g       += dictionaries_to_append[28]
                names_list_light_g       += dictionaries_to_append[29]
                names_list_heavy_h       += dictionaries_to_append[30]
                names_list_light_h       += dictionaries_to_append[31]
                Chem_con                 += dictionaries_to_append[32]

                #print(All_positions_and_chains)
                #ADCs += dictionaries_to_append[10]


##Get labels and positions
            if "H[" not in keyslist[i] and "Linker[" not in keyslist[i] and "X" not in keyslist[i] and "C[" not in keyslist[i]:
                print("checkpoint27")
                Label_Locations = getcoordinates[3]
                text_coordinates.append([Label_Locations])
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [getcoordinates[0], righthanded,direction]
                Domain_Text.append(str(Domain_name)+mod_label)


            elif ("X[" in keyslist[i] or "C[" in keyslist[i]):
                if mod != "Leucine" and "X[" in keyslist[i]:
                    ADCs.append(getcoordinates[0])
                elif mod != "Leucine" and "C[" in keyslist[i]:
                    Chem_con.append(getcoordinates[0])
                if slant == True and righthanded == False:
                    Label_Locations = [getcoordinates[3][0]-20,getcoordinates[3][1]]
                elif slant == True and  righthanded== True:
                    Label_Locations = [getcoordinates[3][0]+20,getcoordinates[3][1]]
                elif slant == False:
                    Label_Locations = [getcoordinates[3][0],getcoordinates[3][1]]
                text_coordinates.append([Label_Locations])
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                Domain_Text.append(str(Domain_name)+mod_label)
                Domain_number = Domain_name+str([i])
                All_positions_and_chains[text] = [getcoordinates[0], righthanded, direction]


##append interactions

            if dictionary != VHa_chain and dictionary != VHb_chain and dictionary != VLa_chain and dictionary != VLb_chain and dictionary != fragment1 and dictionary != fragment2 and dictionary != fragment3 and dictionary != fragment4:

                if len(dictionary.get(keyslist[i])) > 1 and interaction_counter ==0:
                    if dictionary == VHa_1_test:
                        keyslistb = list(VHb_chain.keys())
                        for x in range(len(keyslistb)):
                            if dictionary.get(keyslist[i])[1] == VHb_chain.get(keyslistb[x])[0][0]:
                                if "H[" in str(keyslist[i]) and i+1 == len(keyslist):
                                    if righthanded == False:
                                        first_interaction =[getcoordinates[1][0],getcoordinates[1][1]-40]
                                    elif righthanded == True:
                                        first_interaction =[getcoordinates[1][0],getcoordinates[1][1]-40]
                                else:
                                    first_interaction =[getcoordinates[0][6],getcoordinates[0][7]]
                                    interaction_counter +=1
                    elif dictionary == VHb_1_test:
                        keyslista = list(VHa_chain.keys())
                        for x in range(len(keyslista)):
                            if dictionary.get(keyslist[i])[1] == VHa_chain.get(keyslista[x])[0][0]:
                                if "H[" in str(keyslist[i]) and i+1 == len(keyslist):
                                    if righthanded == False:
                                        first_interaction =[getcoordinates[1][0],getcoordinates[1][1]-40]
                                    elif righthanded == True:
                                        first_interaction =[getcoordinates[1][0],getcoordinates[1][1]-40]
                                else:
                                    first_interaction =[getcoordinates[0][6],getcoordinates[0][7]]
                                    interaction_counter +=1
                    elif dictionary == VLa_1_test:
                        keyslista = list(VHa_1_test.keys())
                        for x in range(len(keyslista)):
                            if dictionary.get(keyslist[i])[1] == VHa_1_test.get(keyslista[x])[0]:
                                first_interaction =[getcoordinates[0][0],getcoordinates[0][1]]
                                interaction_counter +=1
                    elif dictionary == VLb_1_test:
                        keyslistb = list(VHb_1_test.keys())
                        for x in range(len(keyslistb)):
                            if dictionary.get(keyslist[i])[1] == VHb_1_test.get(keyslistb[x])[0]:
                                first_interaction =[getcoordinates[0][0],getcoordinates[0][1]]
                                interaction_counter +=1


                elif len(dictionary.get(keyslist[i])) == 1 and interaction_counter ==0 and ("X[" in keyslist[i] or "C[" in keyslist[i]):
                    if dictionary == VHa_1_test:
                        keyslistb = list(VHb_chain.keys())
                        for x in range(len(keyslistb)):
                            if dictionary.get(keyslist[i])[0] == VHb_chain.get(keyslistb[x])[0][0]:
                                first_interaction =[getcoordinates[0][6],getcoordinates[0][7]]
                                interaction_counter +=1
                    elif dictionary == VHb_1_test:
                        keyslista = list(VHa_chain.keys())
                        for x in range(len(keyslista)):
                            if dictionary.get(keyslist[i])[0] == VHa_chain.get(keyslista[x])[0][0]:
                                first_interaction =[getcoordinates[0][6],getcoordinates[0][7]]
                                interaction_counter +=1
##Protein multimers
            if "X" in keyslist[i] and dictionary.get(keyslist[i])[0] in multimers_keyslist:
                index = dictionary.get(keyslist[i])[0]
                current_multimer = multimers.get(index)[1]
                finished_multimers_numbers.append(current_multimer)
                finished_multimers_indexes.append(index)
                #multimers[dictionary.get(keyslist[i])[0]] = getcoordinates[0]
                finished_multimers.append(getcoordinates[0])


##Sort extra disulphide bridges
            if  dictionary == VHa_chain or dictionary == VHb_chain or dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2 or dictionary == fragment3 or dictionary == fragment4:
                if disulphide_bridge_count > 0 and "H[" not in keyslist[i]:#  and "X" not in keyslist[i]:
                    bottom_bond_extra = getcoordinates[1]
                    bottomx=bottom_bond_extra[0]
                    bottomy=bottom_bond_extra[1]
                    if slant == True and righthanded == False:
                        topx = bottomx+20
                    elif slant == True and righthanded == True:
                        topx = bottomx-20
                    elif slant==False:
                        topx = bottomx
                    topy = getcoordinates[2][1]+40
                    location   = dictionary.get(keyslist[i])[0]
                    interactor = dictionary.get(keyslist[i])[1]
                    H_bond_coordinates = extra_disulphide_maker(disulphide_bridge_count,bottomx,bottomy,topx,topy,righthanded)
                    extra_disulphide_bridges[interactor]=H_bond_coordinates
                    extra_disulphidebridge_keyslist = list(extra_disulphide_bridges.keys())
                    for j in range(len(extra_disulphide_bridges)):
                        if int(extra_disulphidebridge_keyslist[j]) == int(location):
                            for x in range(len(extra_disulphide_bridges.get(interactor))):
                                coordinate1 = extra_disulphide_bridges.get(interactor)[x][0]
                                coordinate2 = extra_disulphide_bridges.get(interactor)[x][1]
                                coordinate3 = extra_disulphide_bridges.get(extra_disulphidebridge_keyslist[j])[x][0]
                                coordinate4 = extra_disulphide_bridges.get(extra_disulphidebridge_keyslist[j])[x][1]
                                coordinates = [coordinate1,coordinate2,coordinate3,coordinate4]
                                completed_disulphidebridges.append(coordinates)

            if "H[" in keyslist[i]:
                bottomx=bottom_bond[0]
                bottomy=bottom_bond[1]
                if H_count %2 !=0:
                    if slant == True and righthanded == False:
                        topx = bottomx+20
                    elif slant == True and righthanded == True:
                        topx = bottomx-20
                    elif slant==False:
                        topx = bottomx
                    topy = bottomy+40
                elif H_count %2 == 0:
                    slant = True
                    if slant == True and righthanded == False:
                        topx = bottomx-20
                    elif slant == True and righthanded == True:
                        topx = bottomx+20
                    elif slant==False:
                        topx = bottomx
                    topy = bottomy+40
                if i+1 != len(dictionary):
                    if Extra_bond == True:
                        Label_Locations = [getcoordinates[3][0],getcoordinates[3][1]+80]
                    else:
                        Label_Locations = [getcoordinates[3][0],getcoordinates[3][1]+60]
                elif i+1 == len(dictionary):
                    if slant == True:
                        top_bond[1] = top_bond[1]+5
                        hinges.append([bottomx,bottomy,topx,topy])
                    elif slant == False:
                        if tangle_found == True and dictionary == VHa_chain:
                            pass
                        else:
                            hinges.append([bottomx,bottomy,topx,topy+20])

                    if Extra_bond == True and righthanded == False:
                        Label_Locations = [getcoordinates[3][0]-60,getcoordinates[3][1]]
                    elif Extra_bond == True and righthanded == True:
                        Label_Locations = [getcoordinates[3][0]+60,getcoordinates[3][1]]
                    elif righthanded== False:
                        Label_Locations = [getcoordinates[3][0]-60,getcoordinates[3][1]-17]
                    elif righthanded== True:
                        Label_Locations = [getcoordinates[3][0]+60,getcoordinates[3][1]-17]
                Hx = Label_Locations[0]
                if righthanded == False and (slant == True or chain_count == 2):
                    Hx -= 0
                elif righthanded == True and (slant == True or chain_count ==2):
                    Hx += 0
                else:
                    Hx = Hx
                Hy = Label_Locations[1]
                Label_Locations = [[Hx,Hy]]
                text_coordinates.append(Label_Locations)
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [Label_Locations, righthanded, innie_or_outie, direction]
                Domain_Text.append(str(Domain_name)+mod_label)

                if H_count == 1:
                    slant=False
                elif H_count == 2:
                    slant=True

            bottom_bond = getcoordinates[1]
            Build_up_downlist.append(Build_up)

        chain    = [x for x in chain if x != []]
        allbonds = bonds
        allbonds = [x for x in bonds if x != []]
        #allbonds = bonds
        return(coordinates_list_heavy_a,coordinates_list_light_a,coordinates_list_heavy_b,coordinates_list_light_b,coordinates_list_heavy_c,coordinates_list_light_c,coordinates_list_heavy_d,coordinates_list_light_d,allbonds,Location_Text,text_coordinates,disulphidebridge1, disulphidebridge2 ,disulphidebridge3,disulphidebridge4,disulphidebridge5,completed_disulphidebridges,Domain_Text,Notes,Notes_positions, arcs_left,arcs_right, arcs_left_slant, arcs_right_slant,ADCs,first_interaction, hinges, linkers,names_list_heavy_a,names_list_light_a,names_list_heavy_b,names_list_light_b,names_list_heavy_c,names_list_light_c,names_list_heavy_d,names_list_light_d,coordinates_list_heavy_e,coordinates_list_light_e,coordinates_list_heavy_f,coordinates_list_light_f,coordinates_list_heavy_g,coordinates_list_light_g,coordinates_list_heavy_h,coordinates_list_light_h,names_list_heavy_e,names_list_light_e,names_list_heavy_f,names_list_light_f,names_list_heavy_g,names_list_light_g,names_list_heavy_h,names_list_light_h,Chem_con)

##Get starting positions
    VHa_2_test = VHa_chain.copy()
    VLa_2_test = VLa_chain.copy()
    VHb_2_test = VHb_chain.copy()
    VLb_2_test = VLb_chain.copy()
    ##tangle detector
    tangle_found = False
    if chain_count ==2:
        keyslista = list(VHa_chain.keys())
        keyslistb = list(VHb_chain.keys())
        try:
            interactor = VHa_chain.get(keyslista[0])[0][1]
            if "Linker[" in keyslista[1]:
                interactor_2 = VHa_chain.get(keyslista[2])[0][1]
            else:
                interactor_2 = VHa_chain.get(keyslista[1])[0][1]
            if interactor_2 < interactor and "H[" in keyslista[-1]:
                tangle_found = True
        except IndexError:
            pass
    ##Get canvas sizes
    width = lower_canvas.winfo_width()
    height = lower_canvas.winfo_height()

    if chain_count == 1:
        VHa_startx, VHa_starty = (width/2),(height/2)-200
        VHb_startx, VHb_starty = 0,0
        VLa_startx, VLa_starty = 0,0
        VLb_startx, VLb_starty = 0,0
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)

    elif  chain_count >= 3 or (chain_count == 2  and tangle_found == False and ("H[" in str(VHa_chain) and "H[" in str(VHb_chain)) or ("X[" in str(VHa_chain) and "X[" in str(VHb_chain)) or ("C[" in str(VHa_chain) and "C[" in str(VHb_chain))):
        VHa_1_test = VHa_chain.copy()
        VLa_1_test = VLa_chain.copy()
        VHb_1_test = VHb_chain.copy()
        VLb_1_test = VLb_chain.copy()



##VHa_chain

        if "H[" in str(VHa_1_test) and "H[" in str(VHb_1_test):
            if IgG2 == False:
                VHa_H_coordinatesx = (width/2)-50
                VHa_H_coordinatesy = (height/2)-50
                VHb_H_coordinatesx = (width/2)+50
                VHb_H_coordinatesy = (height/2)-50
            elif IgG2 == True:
                VHa_H_coordinatesx = (width/2)-200
                VHa_H_coordinatesy = (height/2)-100
                VHb_H_coordinatesx = (width/2)-100
                VHb_H_coordinatesy = (height/2)-100
        elif (("X[" in str(VHa_1_test) and "X[" in str(VHb_1_test)) or ("C[" in str(VHa_1_test) and "C[" in str(VHb_1_test)))and ("H[" not in str(VHa_1_test) and "H[" not in str(VHb_1_test)):
            keyslista = list(VHa_chain.keys())
            keyslistb = list(VHb_chain.keys())
            Xa = ""
            Xb = ""
            for i in range(len(keyslista)):
                if "X[" in keyslista[i] or "C[" in keyslista[i] :
                    Xa = int(VHa_chain.get(keyslista[i])[0][0])
                    for i in range(len(keyslistb)):
                        if "X[" in keyslistb[i] or "C[" in keyslistb[i]:
                            Xb = int(VHb_chain.get(keyslistb[i])[0][0])
                            if Xa == Xb:
                                VHa_H_coordinatesx = (width/2)
                                VHa_H_coordinatesy = (height/2)-50
                                VHb_H_coordinatesx = (width/2)
                                VHb_H_coordinatesy = (height/2)-50
                            elif Xa != Xb:
                                VHa_H_coordinatesx = (width/2)-50
                                VHa_H_coordinatesy = (height/2)-50
                                VHb_H_coordinatesx = (width/2)+50
                                VHb_H_coordinatesy = (height/2)-50

        elif chain_count == 4 and "H[" not in str(VHa_1_test) and "H[" not in str(VHb_1_test) and "X[" not in str(VHa_1_test) and "X[" not in str(VHb_1_test) and "C[" not in str(VHa_1_test) and "C[" not in str(VHb_1_test):
            VHa_H_coordinatesx = (width/2)-32
            VHa_H_coordinatesy = (height/2)
            VHb_H_coordinatesx = (width/2)+32
            VHb_H_coordinatesy = (height/2)

        teststartx = 0
        teststarty = 0
        testHpositionVHa = renderchains(VHa_1_test,teststartx,teststarty)[25]
        test_H_positionx = testHpositionVHa[0]
        test_H_positiony = testHpositionVHa[1]
        differencetest_desiredx = VHa_H_coordinatesx - test_H_positionx
        differencetest_desiredy = VHa_H_coordinatesy - test_H_positiony
        VHa_startx = teststartx + differencetest_desiredx
        VHa_starty = teststarty + differencetest_desiredy

##VHb_chain

        teststartx = 800
        teststarty = 0
        testHpositionVHb = renderchains(VHb_1_test,teststartx,teststarty)[25]
        test_H_positionx = testHpositionVHb[0]
        test_H_positiony = testHpositionVHb[1]
        differencetest_desiredx = test_H_positionx - VHb_H_coordinatesx
        differencetest_desiredy = VHb_H_coordinatesy - test_H_positiony
        VHb_startx = teststartx - differencetest_desiredx
        VHb_starty = teststarty + differencetest_desiredy


        All_positions_and_chains    ={}
        extra_disulphide_bridges    ={}
        H_disulphide_coordinates    ={}
        completed_disulphidebridges=[]
        Notes           = []
        Notes_positions = []
        noted_Hbonds=False
        H_coordinatey   = []
        #print("MULTIMER RESET")
        multimers =  protein_multimer(VHa_chain,VLa_chain,VHb_chain,VLb_chain,fragment1,fragment2,fragment3,fragment4)
        multimers_keyslist = (multimers.keys())
        finished_multimers = []
        finished_multimers_numbers = []
        finidhed_multimers_indexes = []



###Render Heavy chains
        if VHb_startx > width and IgG2 == True:
            VHb_startx -= 515
        elif VHb_startx > width and IgG2 == False:
            VHb_startx -= 600
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)


###Get start positions of light chains and render
    elif chain_count == 2 :

        keyslista = list(VHa_chain.keys())
        keyslistb = list(VHb_chain.keys())
        keyslist = list(VHa_chain.keys())
        VHb_startx, VHb_starty = (width/2)+100,(height/2)-200
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VHa_list = list(VHa_chain.keys())
        VHa_startx, VHa_starty= (width/2)-150,(height/2)-200
        try:
            VHa_inter = VHa_chain.get(VHa_list[0])[0][1]
            if VHa_chain.get(keyslista[0])[0][0] == VHb_chain.get(keyslistb[0])[1] or VHa_chain.get(keyslista[0])[0][1] in All_positions_and_chains:
                VHa_start = find_the_fragment(VHa_inter,All_positions_and_chains)
                if VHa_start is not None:
                    righthanded = VHa_start[1]
                    if righthanded == True:
                        VHa_startx=VHa_start[0][0]-60
                        VHa_starty=VHa_start[0][1]
                    elif righthanded==False:
                        VHa_startx=VHa_start[0][0]-60
                        VHa_starty=VHa_start[0][1]
                else:
                    VHa_startx, VHa_starty= (width/2)-150,(height/2)-200
        except IndexError:
            if len(VHa_chain.get(VHa_list[1])[0]) > 1:
                if VHa_chain.get(VHa_list[1])[0][0] == VHb_chain.get(keyslistb[1])[1]:
                    VHa_startx, VHa_starty= (width/2)-50,(height/2)-200

        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)



##Light chains_list

    if VLa_chain != {}:
        All_positions_and_chains_list = list(All_positions_and_chains.keys())
        keyslistlight_a = list(VLa_chain.keys())
        interactor_count = 0
        if interactor_count == 0:
            for i in range(len(keyslistlight_a)):
                if len(VLa_chain.get(keyslistlight_a[i])[0]) >1:
                    interactor = VLa_chain.get(keyslistlight_a[i])[0][1]
                else:
                    continue
                for j in range(len(All_positions_and_chains_list)):
                    if interactor == All_positions_and_chains_list[j] and interactor_count == 0:
                        interactor_count +=1
                        VLa_reference_pos = find_the_fragment2(interactor,All_positions_and_chains)
                        VLa_refx = VLa_reference_pos[0]
                        VLa_refy = VLa_reference_pos[1]
                        teststartx = 0
                        teststarty = 0
                        testHpositionVLa = renderchains(VLa_1_test,teststartx,teststarty)[25]
                        test_H_positionx = testHpositionVLa[0]
                        test_H_positiony = testHpositionVLa[1]
                        differencetest_desiredx = VLa_refx - test_H_positionx
                        differencetest_desiredy = VLa_refy - test_H_positiony
                        VLa_startx = teststartx + differencetest_desiredx
                        VLa_starty = teststarty + differencetest_desiredy
                        break
    else:
        VLa_startx, VLa_starty = 0,0

    if VLb_chain != {}:
        All_positions_and_chains_list = list(All_positions_and_chains.keys())
        keyslistlight_b = list(VLb_chain.keys())
        interactor_count = 0
        if interactor_count == 0:
            for i in range(len(keyslistlight_b)):
                if len(VLb_chain.get(keyslistlight_b[i])[0]) >1:
                    interactor = VLb_chain.get(keyslistlight_b[i])[0][1]
                else:
                    continue
                for j in range(len(All_positions_and_chains_list)):
                    if interactor == All_positions_and_chains_list[j] and interactor_count == 0:
                        interactor_count += 1
                        VLb_reference_pos = find_the_fragment2(interactor,All_positions_and_chains)
                        VLb_refx = VLb_reference_pos[0]
                        VLb_refy = VLb_reference_pos[1]
                        teststartx = 800
                        teststarty = 0
                        testHpositionVLb = renderchains(VLb_1_test,teststartx,teststarty)[25]
                        test_H_positionx = testHpositionVLb[0]
                        test_H_positiony = testHpositionVLb[1]
                        differencetest_desiredx = VLb_refx - test_H_positionx
                        differencetest_desiredy = VLb_refy - test_H_positiony
                        VLb_startx = teststartx + differencetest_desiredx
                        VLb_starty = teststarty + differencetest_desiredy
                        break
    else:
        VLb_startx,VLb_starty = 0,0

    VLa_stats = renderchains(VLa_chain,VLa_startx, VLa_starty)
    VLb_stats = renderchains(VLb_chain,VLb_startx, VLb_starty)


##Two-chained Abs


    frag1_startx,frag1_starty,frag2_startx,frag2_starty,frag3_startx,frag3_starty,frag4_startx,frag4_starty = 0,0,0,0,0,0,0,0
    if IgG2 == False:
        if fragment1 != {}:
            fragment1_list = list(fragment1.keys())
            fragment_inter = fragment1.get(fragment1_list[0])[0][1]
            frag1_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag1_start[1]
            if righthanded == True:
                frag1_startx=frag1_start[0][0]+60
                frag1_starty=frag1_start[0][1]
            elif righthanded==False:
                frag1_startx=frag1_start[0][0]-60
                frag1_starty=frag1_start[0][1]
        if fragment2 != {}:
            fragment2_list = list(fragment2.keys())
            fragment_inter = fragment2.get(fragment2_list[0])[0][1]
            frag2_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag2_start[1]
            if righthanded == True:
                frag2_startx=frag2_start[0][0]+60
                frag2_starty=frag2_start[0][1]
            elif righthanded==False:
                frag2_startx=frag2_start[0][0]-60
                frag2_starty=frag2_start[0][1]
        if fragment3 != {}:
            fragment3_list = list(fragment3.keys())
            fragment_inter = fragment3.get(fragment3_list[0])[0][1]
            frag3_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag3_start[1]
            if righthanded == True:
                frag3_startx=frag3_start[0][0]+60
                frag3_starty=frag3_start[0][1]
            elif righthanded==False:
                frag3_startx=frag3_start[0][0]-60
                frag3_starty=frag3_start[0][1]
        if fragment4 != {}:
            fragment4_list = list(fragment4.keys())
            fragment_inter = fragment4.get(fragment4_list[0])[0][1]
            frag4_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag4_start[1]
            if righthanded == True:
                frag4_startx=frag4_start[0][0]+60
                frag4_starty=frag4_start[0][1]
            elif righthanded==False:
                frag4_startx=frag4_start[0][0]-60
                frag4_starty=frag4_start[0][1]
        frag1_stat= renderchains(fragment1,frag1_startx,frag1_starty)
        frag2_stat= renderchains(fragment2,frag2_startx,frag2_starty)
        frag3_stat= renderchains(fragment3,frag3_startx,frag3_starty)
        frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty)
    elif IgG2 == True:
        if (("X" in str(list(fragment1.keys())) and "X" in str(list(fragment3.keys()))) or ("C[" in str(list(fragment1.keys())) and "C[" in str(list(fragment3.keys())))) and ("CH2" not in str(list(fragment1.keys())) and "CH2" not in str(list(fragment3.keys()))) :
            frag1_stat= renderchains(fragment1,VHa_startx-0,VHa_starty+200)
            test_H_positionVHa = frag1_stat[25]
            test_H_positionx = testHpositionVHa[0]
            test_H_positiony = testHpositionVHa[1]
            differencetest_desiredx = VHa_H_coordinatesx - test_H_positionx
            differencetest_desiredy = VHa_H_coordinatesy - test_H_positiony
            #
            #
            fragment2_list = list(fragment2.keys())
            fragment_inter = fragment2.get(fragment2_list[0])[0][1]
            frag2_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag2_start[1]
            if righthanded == True:
                frag2_startx=frag2_start[0][0]+60
                frag2_starty=frag2_start[0][1]
            elif righthanded==False:
                frag2_startx=frag2_start[0][0]-60
                frag2_starty=frag2_start[0][1]
            frag2_stat = renderchains(fragment2,frag2_startx,frag2_starty)

            frag3_stat= renderchains(fragment3,VHb_startx+0,VHb_starty+200)
            fragment4_list = list(fragment4.keys())
            fragment_inter = fragment4.get(fragment4_list[0])[0][1]
            frag4_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag4_start[1]
            if righthanded == True:
                frag4_startx=frag4_start[0][0]+60
                frag4_starty=frag4_start[0][1]
            elif righthanded==False:
                frag4_startx=frag4_start[0][0]-60
                frag4_starty=frag4_start[0][1]
            frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty)

        else:

            VHa_H_coordinatesx = ((width/6)*4)-10
            VHa_H_coordinatesy = (height/2)+100
            VHb_H_coordinatesx = VHa_H_coordinatesx+100
            VHb_H_coordinatesy = (height/2)+100

            frag1_stat= renderchains(fragment1,VHa_startx+325,VHa_starty+200)
            test_H_positionfrag1 = frag1_stat[26]


            test_H_positionx = test_H_positionfrag1[0][0]
            test_H_positiony = test_H_positionfrag1[0][1]
            frag1_differencetest_desiredx = test_H_positionx - VHa_H_coordinatesx
            frag1_differencetest_desiredy = test_H_positiony - VHa_H_coordinatesy
            coordinates_to_change = [frag1_stat[0],frag1_stat[1],frag1_stat[2],frag1_stat[3],frag1_stat[4],frag1_stat[5],frag1_stat[6],frag1_stat[7],frag1_stat[36],frag1_stat[37],frag1_stat[38],frag1_stat[39],frag1_stat[40],frag1_stat[41],frag1_stat[42],frag1_stat[43],frag1_stat[8],frag1_stat[26],frag1_stat[27],frag1_stat[24],frag1_stat[11], frag1_stat[10]]
            conj_x1 = get_min_max_coordinates(frag1_stat[52][0])[0]
            conj_x2 = get_min_max_coordinates(frag1_stat[52][0])[1]
            conj_y1 = get_min_max_coordinates(frag1_stat[52][0])[2]
            conj_y2 = get_min_max_coordinates(frag1_stat[52][0])[3]
            conj_fixed = False
            for i in range(len(coordinates_to_change)):
                for j in range(len(coordinates_to_change[i])):
                    for k in range(len(coordinates_to_change[i][j])):
                        if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
                            try:
                                if conj_x1 <= coordinates_to_change[i][j][2] <= conj_x2 and conj_y1 <= coordinates_to_change[i][j][3] <= conj_y2 and conj_fixed == False:
                                    print("LINK TO CONJUGATE")
                                    coordinates_to_change[i][j][0]-= frag1_differencetest_desiredx
                                    coordinates_to_change[i][j][1]+= frag1_differencetest_desiredy
                                    coordinates_to_change[i][j][2]= frag1_stat[52][0][0] + frag1_differencetest_desiredx
                                    coordinates_to_change[i][j][3]= frag1_stat[52][0][1] + frag1_differencetest_desiredy+20
                                    conj_fixed = True

#
                                else:
                                    if k %2 != 0:
                                        coordinates_to_change[i][j][k] -= frag1_differencetest_desiredy
                                    elif k%2 ==0:
                                        coordinates_to_change[i][j][k] -= frag1_differencetest_desiredx

                            except IndexError:
                                if  k %2 != 0:
                                    coordinates_to_change[i][j][k] -= frag1_differencetest_desiredy
                                elif k %2 == 0:
                                    coordinates_to_change[i][j][k] -= frag1_differencetest_desiredx
                        else:
                            for l in range(len(coordinates_to_change[i][j][k])):
                                if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
                                    if l %2 != 0:
                                        coordinates_to_change[i][j][k][l] -= frag1_differencetest_desiredy
                                    elif l%2 ==0:
                                        coordinates_to_change[i][j][k][l] -= frag1_differencetest_desiredx
            #All_positions_and_chains_list = list(All_positions_and_chains.keys())
            #print(All_positions_and_chains)
            for i in range(len(frag1_stat[52][0])):
                if i % 2 == 0:
                    frag1_stat[52][0][i] += frag1_differencetest_desiredx
                elif i % 2 != 0:
                    frag1_stat[52][0][i] += frag1_differencetest_desiredy

            fragment1_keyslist = list(fragment1.keys())
            for i in range(len(All_positions_and_chains_list)):
                number = All_positions_and_chains_list[i]
                for j in range(len(fragment1_keyslist)):
                    try:
                        keyslist_number = fragment1.get(fragment1_keyslist[j])[0]
                        if number == keyslist_number:
                            for l in range(len(All_positions_and_chains.get(All_positions_and_chains_list[i])[0])):
                                #if All_positions_and_chains.get(All_positions_and_chains_list[i])[0] != frag1_stat[52][0]:
                                if l %2 != 0:
                                    All_positions_and_chains.get(All_positions_and_chains_list[i])[0][l] -= frag1_differencetest_desiredy
                                elif l%2 ==0:
                                    All_positions_and_chains.get(All_positions_and_chains_list[i])[0][l] -= frag1_differencetest_desiredx
                    except IndexError:
                        pass
            disulphide_keyslist = list(extra_disulphide_bridges.keys())
            for i in range(len(disulphide_keyslist)):
                number = disulphide_keyslist[i]
                for j in range(len(fragment1_keyslist)):
                    try:
                        keyslist_number = fragment1.get(fragment1_keyslist[j])[1]

                        if number == keyslist_number:

                            for l in range(len(extra_disulphide_bridges.get(disulphide_keyslist[i])[0])):
                                if l %2 != 0:
                                    extra_disulphide_bridges.get(disulphide_keyslist[i])[0][l] -= frag1_differencetest_desiredy
                                elif l%2 ==0:
                                    extra_disulphide_bridges.get(disulphide_keyslist[i])[0][l] -= frag1_differencetest_desiredx
                    except IndexError:
                        pass
            for i in range(len(frag1_stat[17])):
                if frag1_stat[17][i] == "C":
                    #for j in range(len(frag1_stat[9][i])):
                        #if j % 2 == 0:
                    frag1_stat[10][i][0][0] += frag1_differencetest_desiredx
                        #elif j % 2 != 0:
                    frag1_stat[10][i][0][1] += frag1_differencetest_desiredy


            #H_keyslist = list(H_disulphide_coordinates.keys())
            #for i in range(len(H_keyslist)):
            #    number = H_keyslist[i]
            #    for j in range(len(fragment1_keyslist)):
            #        try:
            #            keyslist_number = fragment1.get(fragment1_keyslist[j])[1]
            #            print()
            #            if number == keyslist_number:
            #                for l in range(len(H_disulphide_coordinates.get(H_keyslist[i])[0])):
            #                    print("ELLO ELLO", H_disulphide_coordinates.get(H_keyslist[i])[0][l])
            #                    if l %2 != 0:
            #                        H_disulphide_coordinates.get(H_keyslist[i])[0][l] += frag1_differencetest_desiredy
            #                    elif l%2 ==0:
            #                        H_disulphide_coordinates.get(H_keyslist[i])[0][l] += frag1_differencetest_desiredx
            #        except IndexError:
            #            pass
            fragment2_list = list(fragment2.keys())
            fragment_inter = fragment2.get(fragment2_list[0])[0][1]
            frag2_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag2_start[1]
            if righthanded == True:
                frag2_startx=frag2_start[0][0]+60
                frag2_starty=frag2_start[0][1]
            elif righthanded==False:
                frag2_startx=frag2_start[0][0]-60
                frag2_starty=frag2_start[0][1]
            frag2_stat = renderchains(fragment2,frag2_startx,frag2_starty)

            frag3_stat= renderchains(fragment3,VHb_startx+325,VHb_starty+200)
            test_H_positionfrag3 = frag3_stat[26]


            test_H_positionx = test_H_positionfrag3[0][0]
            test_H_positiony = test_H_positionfrag3[0][1]
            frag3_differencetest_desiredx = VHb_H_coordinatesx- test_H_positionx
            frag3_differencetest_desiredy = VHb_H_coordinatesy- test_H_positiony
            coordinates_to_change = [frag3_stat[0],frag3_stat[1],frag3_stat[2],frag3_stat[3],frag3_stat[4],frag3_stat[5],frag3_stat[6],frag3_stat[7],frag3_stat[36],frag3_stat[37],frag3_stat[38],frag3_stat[39],frag3_stat[40],frag3_stat[41],frag3_stat[42],frag3_stat[43],frag3_stat[8],frag3_stat[26],frag3_stat[27],frag3_stat[24],frag3_stat[11], frag3_stat[10]]
            for i in range(len(coordinates_to_change)):
                for j in range(len(coordinates_to_change[i])):
                    for k in range(len(coordinates_to_change[i][j])):
                        print(coordinates_to_change[i][j][k])
                        if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:

                            if  k %2 != 0:
                                coordinates_to_change[i][j][k] += frag3_differencetest_desiredy
                            elif k %2 == 0:
                                coordinates_to_change[i][j][k] += frag3_differencetest_desiredx
                        else:
                            for l in range(len(coordinates_to_change[i][j][k])):
                                if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
                                    if l %2 != 0:
                                        coordinates_to_change[i][j][k][l] += frag3_differencetest_desiredy
                                    elif l%2 ==0:
                                        coordinates_to_change[i][j][k][l] += frag3_differencetest_desiredx
            fragment3_keyslist = list(fragment3.keys())
            disulphide_keyslist = list(extra_disulphide_bridges.keys())
            for i in range(len(disulphide_keyslist)):
                number = disulphide_keyslist[i]
                for j in range(len(fragment3_keyslist)):
                    try:
                        keyslist_number = fragment3.get(fragment3_keyslist[j])[1]
                        if number == keyslist_number:

                            for l in range(len(extra_disulphide_bridges.get(disulphide_keyslist[i])[0])):
                                if l %2 != 0:
                                    extra_disulphide_bridges.get(disulphide_keyslist[i])[0][l] += frag3_differencetest_desiredy
                                elif l%2 ==0:
                                    extra_disulphide_bridges.get(disulphide_keyslist[i])[0][l] += frag3_differencetest_desiredx
                    except IndexError:
                        pass



            #H_keyslist = list(H_disulphide_coordinates.keys())
            #for i in range(len(H_keyslist)):
            #    number = H_keyslist[i]
            #    for j in range(len(fragment3_keyslist)):
            #        try:
            #            keyslist_number = fragment3.get(fragment3_keyslist[j])[1]
            #            if number == keyslist_number:
            #                for l in range(len(H_disulphide_coordinates.get(H_keyslist[i])[0])):
            #                    print("ELLO ELLO", H_disulphide_coordinates.get(H_keyslist[i])[0][l])
            #                    if l %2 != 0:
            #                        H_disulphide_coordinates.get(H_keyslist[i])[0][l] += differencetest_desiredy
            #                    elif l%2 ==0:
            #                        H_disulphide_coordinates.get(H_keyslist[i])[0][l] += differencetest_desiredx
            #        except IndexError:
            #            pass

            fragment4_list = list(fragment4.keys())
            fragment_inter = fragment4.get(fragment4_list[0])[0][1]
            frag4_start = find_the_fragment(fragment_inter,All_positions_and_chains)
            righthanded = frag4_start[1]
            if righthanded == True:
                frag4_startx=frag4_start[0][0]+60
                frag4_starty=frag4_start[0][1]
            elif righthanded==False:
                frag4_startx=frag4_start[0][0]-60
                frag4_starty=frag4_start[0][1]
            frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty)
            hingey1 = frag3_stat[26][0][1]


            for i in range(len(completed_disulphidebridges)):
                if completed_disulphidebridges[i][1] > hingey1:
                    completed_disulphidebridges[i][0] -= frag1_differencetest_desiredx
                    completed_disulphidebridges[i][1] -= frag1_differencetest_desiredy
                    completed_disulphidebridges[i][2] += frag3_differencetest_desiredx
                    completed_disulphidebridges[i][3] += frag3_differencetest_desiredy


    keyslist_All_positions_and_chains = list(All_positions_and_chains.keys())
    Yout_of_range = False
    Yhow_much = 0
    Xout_of_range = False
    Xhow_much = 0
    for i in range(len(keyslist_All_positions_and_chains)):
        for j in range(len(All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0])):
            if isinstance(All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j], float) == True:
                if All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j] < 0:
                    Yout_of_range = True
                    if Yhow_much > All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j]:
                        Yhow_much = All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j]
    for i in range(len(keyslist_All_positions_and_chains)):
        for j in range(len(All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0])):
            if isinstance(All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j], float) == True:
                if All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j] > 800:
                    Xout_of_range = True
                    if Xhow_much < All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j]:
                        Xhow_much = All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j]

    if Yout_of_range == True:
        Ynew_start = Yhow_much+5
    #if Xout_of_range == True:
    #    Xnew_start = Xhow_much-1000
    #print("X out og range", Xnew_start)




    Heavy_Domains_a     = VHa_stats[0] + VLa_stats[0] + VHb_stats[0] + VLb_stats[0] + frag1_stat[0] + frag2_stat[0] + frag3_stat[0] + frag4_stat[0]
    Light_Domains_a     = VHa_stats[1] + VLa_stats[1] + VHb_stats[1] + VLb_stats[1] + frag1_stat[1] + frag2_stat[1] + frag3_stat[1] + frag4_stat[1]
    Heavy_Domains_b     = VHa_stats[2] + VLa_stats[2] + VHb_stats[2] + VLb_stats[2] + frag1_stat[2] + frag2_stat[2] + frag3_stat[2] + frag4_stat[2]
    Light_Domains_b     = VHa_stats[3] + VLa_stats[3] + VHb_stats[3] + VLb_stats[3] + frag1_stat[3] + frag2_stat[3] + frag3_stat[3] + frag4_stat[3]
    Heavy_Domains_c     = VHa_stats[4] + VLa_stats[4] + VHb_stats[4] + VLb_stats[4] + frag1_stat[4] + frag2_stat[4] + frag3_stat[4] + frag4_stat[4]
    Light_Domains_c     = VHa_stats[5] + VLa_stats[5] + VHb_stats[5] + VLb_stats[5] + frag1_stat[5] + frag2_stat[5] + frag3_stat[5] + frag4_stat[5]
    Heavy_Domains_d     = VHa_stats[6] + VLa_stats[6] + VHb_stats[6] + VLb_stats[6] + frag1_stat[6] + frag2_stat[6] + frag3_stat[6] + frag4_stat[6]
    Light_Domains_d     = VHa_stats[7] + VLa_stats[7] + VHb_stats[7] + VLb_stats[7] + frag1_stat[7] + frag2_stat[7] + frag3_stat[7] + frag4_stat[7]
    Bonds               = VHa_stats[8] + VLa_stats[8] + VHb_stats[8] + VLb_stats[8] + frag1_stat[8] + frag2_stat[8] + frag3_stat[8] + frag4_stat[8]
    Bonds               = [x for x in Bonds if x != []]
    disulphide_bridges  = [VHa_stats[11] + VHb_stats[11]] + [VHa_stats[12] + VHb_stats[12]] + [VHa_stats[13] + VHb_stats[13]] + [VHa_stats[14] + VHb_stats[14]]+ [VHa_stats[15] + VHb_stats[15]] + completed_disulphidebridges
    disulphide_bridges  = [x for x in disulphide_bridges if x != []]
    Label_Text          = VHa_stats[9] + VLa_stats[9] + VHb_stats[9] + VLb_stats[9] + frag1_stat[9] + frag2_stat[9] + frag3_stat[9] + frag4_stat[9]
    Label_spot          = VHa_stats[10] + VLa_stats[10] + VHb_stats[10] + VLb_stats[10] + frag1_stat[10] + frag2_stat[10] + frag3_stat[10] + frag4_stat[10]
    Hinges              = VHa_stats[26] + VLa_stats[26] + VHb_stats[26] + VLb_stats[26] + frag1_stat[26] + frag2_stat[26] + frag3_stat[26] + frag4_stat[26]
    Hinges              = [x for x in Hinges if x != []]
    Linkers             = VHa_stats[27] + VLa_stats[27] + VHb_stats[27] + VLb_stats[27] + frag1_stat[27] + frag2_stat[27] + frag3_stat[27] + frag4_stat[27]
    Linkers             = [x for x in Linkers if x != []]
    Domain_Text         = VHa_stats[17] + VLa_stats[17] + VHb_stats[17] + VLb_stats[17] + frag1_stat[17] + frag2_stat[17] + frag3_stat[17] + frag4_stat[17]
    Notes               = VHa_stats[18] + VLa_stats[18] + VHb_stats[18] + VLb_stats[18] + frag1_stat[18] + frag2_stat[18] + frag3_stat[18] + frag4_stat[18]
    Notes_positions     = VHa_stats[19] + VLa_stats[19] + VHb_stats[19] + VLb_stats[19] + frag1_stat[19] + frag2_stat[19] + frag3_stat[19] + frag4_stat[19]
    arcs_left           = VHa_stats[20] + VLa_stats[20] + VHb_stats[20] + VLb_stats[20] + frag1_stat[20] + frag2_stat[20] + frag3_stat[20] + frag4_stat[20]
    arcs_right          = VHa_stats[21] + VLa_stats[21] + VHb_stats[21] + VLb_stats[21] + frag1_stat[21] + frag2_stat[21] + frag3_stat[21] + frag4_stat[21]
    arcs_left_slant     = VHa_stats[22] + VLa_stats[22] + VHb_stats[22] + VLb_stats[22] + frag1_stat[22] + frag2_stat[22] + frag3_stat[22] + frag4_stat[22]
    arcs_right_slant    = VHa_stats[23] + VLa_stats[23] + VHb_stats[23] + VLb_stats[23] + frag1_stat[23] + frag2_stat[23] + frag3_stat[23] + frag4_stat[23]
    ADCs                = VHa_stats[24] + VLa_stats[24] + VHb_stats[24] + VLb_stats[24] + frag1_stat[24] + frag2_stat[24] + frag3_stat[24] + frag4_stat[24]
    names_list_heavy_a  = VHa_stats[28] + VLa_stats[28] + VHb_stats[28] + VLb_stats[28] + frag1_stat[28] + frag2_stat[28] + frag3_stat[28] + frag4_stat[28]
    names_list_light_a  = VHa_stats[29] + VLa_stats[29] + VHb_stats[29] + VLb_stats[29] + frag1_stat[29] + frag2_stat[29] + frag3_stat[29] + frag4_stat[29]
    names_list_heavy_b  = VHa_stats[30] + VLa_stats[30] + VHb_stats[30] + VLb_stats[30] + frag1_stat[30] + frag2_stat[30] + frag3_stat[30] + frag4_stat[30]
    names_list_light_b  = VHa_stats[31] + VLa_stats[31] + VHb_stats[31] + VLb_stats[31] + frag1_stat[31] + frag2_stat[31] + frag3_stat[31] + frag4_stat[31]
    names_list_heavy_c  = VHa_stats[32] + VLa_stats[32] + VHb_stats[32] + VLb_stats[32] + frag1_stat[32] + frag2_stat[32] + frag3_stat[32] + frag4_stat[32]
    names_list_light_c  = VHa_stats[33] + VLa_stats[33] + VHb_stats[33] + VLb_stats[33] + frag1_stat[33] + frag2_stat[33] + frag3_stat[33] + frag4_stat[33]
    names_list_heavy_d  = VHa_stats[34] + VLa_stats[34] + VHb_stats[34] + VLb_stats[34] + frag1_stat[34] + frag2_stat[34] + frag3_stat[34] + frag4_stat[34]
    names_list_light_d  = VHa_stats[35] + VLa_stats[35] + VHb_stats[35] + VLb_stats[35] + frag1_stat[35] + frag2_stat[35] + frag3_stat[35] + frag4_stat[35]
    Heavy_Domains_e     = VHa_stats[36] + VLa_stats[36] + VHb_stats[36] + VLb_stats[36] + frag1_stat[36] + frag2_stat[36] + frag3_stat[36] + frag4_stat[36]
    Light_Domains_e     = VHa_stats[37] + VLa_stats[37] + VHb_stats[37] + VLb_stats[37] + frag1_stat[37] + frag2_stat[37] + frag3_stat[37] + frag4_stat[37]
    Heavy_Domains_f     = VHa_stats[38] + VLa_stats[38] + VHb_stats[38] + VLb_stats[38] + frag1_stat[38] + frag2_stat[38] + frag3_stat[38] + frag4_stat[38]
    Light_Domains_f     = VHa_stats[39] + VLa_stats[39] + VHb_stats[39] + VLb_stats[39] + frag1_stat[39] + frag2_stat[39] + frag3_stat[39] + frag4_stat[39]
    Heavy_Domains_g     = VHa_stats[40] + VLa_stats[40] + VHb_stats[40] + VLb_stats[40] + frag1_stat[40] + frag2_stat[40] + frag3_stat[40] + frag4_stat[40]
    Light_Domains_g     = VHa_stats[41] + VLa_stats[41] + VHb_stats[41] + VLb_stats[41] + frag1_stat[41] + frag2_stat[41] + frag3_stat[41] + frag4_stat[41]
    Heavy_Domains_h     = VHa_stats[42] + VLa_stats[42] + VHb_stats[42] + VLb_stats[42] + frag1_stat[42] + frag2_stat[42] + frag3_stat[42] + frag4_stat[42]
    Light_Domains_h     = VHa_stats[43] + VLa_stats[43] + VHb_stats[43] + VLb_stats[43] + frag1_stat[43] + frag2_stat[43] + frag3_stat[43] + frag4_stat[43]
    names_list_heavy_e  = VHa_stats[44] + VLa_stats[44] + VHb_stats[44] + VLb_stats[44] + frag1_stat[44] + frag2_stat[44] + frag3_stat[44] + frag4_stat[44]
    names_list_light_e  = VHa_stats[45] + VLa_stats[45] + VHb_stats[45] + VLb_stats[45] + frag1_stat[45] + frag2_stat[45] + frag3_stat[45] + frag4_stat[45]
    names_list_heavy_f  = VHa_stats[46] + VLa_stats[46] + VHb_stats[46] + VLb_stats[46] + frag1_stat[46] + frag2_stat[46] + frag3_stat[46] + frag4_stat[46]
    names_list_light_f  = VHa_stats[47] + VLa_stats[47] + VHb_stats[47] + VLb_stats[47] + frag1_stat[47] + frag2_stat[47] + frag3_stat[47] + frag4_stat[47]
    names_list_heavy_g  = VHa_stats[48] + VLa_stats[48] + VHb_stats[48] + VLb_stats[48] + frag1_stat[48] + frag2_stat[48] + frag3_stat[48] + frag4_stat[48]
    names_list_light_g  = VHa_stats[49] + VLa_stats[49] + VHb_stats[49] + VLb_stats[49] + frag1_stat[49] + frag2_stat[49] + frag3_stat[49] + frag4_stat[49]
    names_list_heavy_h  = VHa_stats[50] + VLa_stats[50] + VHb_stats[50] + VLb_stats[50] + frag1_stat[50] + frag2_stat[50] + frag3_stat[50] + frag4_stat[50]
    names_list_light_h  = VHa_stats[51] + VLa_stats[51] + VHb_stats[51] + VLb_stats[51] + frag1_stat[51] + frag2_stat[51] + frag3_stat[51] + frag4_stat[51]
    Chem_con            = VHa_stats[52] + VLa_stats[52] + VHb_stats[52] + VLb_stats[52] + frag1_stat[52] + frag2_stat[52] + frag3_stat[52] + frag4_stat[52]
    if Yout_of_range == True:
        print("OH MY GOSH YES")
        Ynew_start = Yhow_much-10
        coordinates_to_change = [Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Heavy_Domains_e,Light_Domains_e,Heavy_Domains_f,Light_Domains_f,Heavy_Domains_g,Light_Domains_g,Heavy_Domains_h,Light_Domains_h,Bonds,Hinges,Linkers,ADCs,disulphide_bridges, Label_spot]
        for i in range(len(coordinates_to_change)):
            for j in range(len(coordinates_to_change[i])):
                for k in range(len(coordinates_to_change[i][j])):
                    print(coordinates_to_change[i][j][k])
                    if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
                        if  k %2 != 0:
                            coordinates_to_change[i][j][k] -= Ynew_start
                    else:
                        for l in range(len(coordinates_to_change[i][j][k])):
                            if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
                                if l %2 != 0:
                                    coordinates_to_change[i][j][k][l] -= Ynew_start
    #if Xout_of_range == True:
    #    Xnew_start = Xhow_much-10
    #    coordinates_to_change = [Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Heavy_Domains_e,Light_Domains_e,Heavy_Domains_f,Light_Domains_f,Heavy_Domains_g,Light_Domains_g,Heavy_Domains_h,Light_Domains_h,Bonds,Hinges,Linkers,ADCs,disulphide_bridges, Label_spot]
    #    for i in range(len(coordinates_to_change)):
    #        for j in range(len(coordinates_to_change[i])):
    #            for k in range(len(coordinates_to_change[i][j])):
    #                print(coordinates_to_change[i][j][k])
    #                if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
    #                    if  k %2 == 0:
    #                        coordinates_to_change[i][j][k] += Xnew_start
    #                else:
    #                    for l in range(len(coordinates_to_change[i][j][k])):
    #                        if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
    #                            if l %2 == 0:
    #                                coordinates_to_change[i][j][k][l] += Xnew_start
    print(Heavy_Domains_a,names_list_heavy_a,Light_Domains_a,names_list_light_a,Heavy_Domains_b,names_list_heavy_b,Light_Domains_b,names_list_light_b,Heavy_Domains_c,names_list_heavy_c,Light_Domains_c,names_list_light_c,Heavy_Domains_d,names_list_heavy_d,Light_Domains_d,names_list_light_d,Label_Text,Label_spot,Domain_Text,Notes,Notes_positions,arcs_left,arcs_right,arcs_left_slant,arcs_right_slant,ADCs,Heavy_Domains_e,names_list_heavy_e,Light_Domains_e,names_list_light_e,Heavy_Domains_f,names_list_heavy_f,Light_Domains_f,names_list_light_f,Heavy_Domains_g,names_list_heavy_g,Light_Domains_g,names_list_light_g,Heavy_Domains_h,names_list_heavy_h,Light_Domains_h,names_list_light_h)
    return(Bonds,disulphide_bridges,Hinges,Linkers,Heavy_Domains_a,names_list_heavy_a,Light_Domains_a,names_list_light_a,Heavy_Domains_b,names_list_heavy_b,Light_Domains_b,names_list_light_b,Heavy_Domains_c,names_list_heavy_c,Light_Domains_c,names_list_light_c,Heavy_Domains_d,names_list_heavy_d,Light_Domains_d,names_list_light_d,Label_Text,Label_spot,Domain_Text,Notes,Notes_positions,arcs_left,arcs_right,arcs_left_slant,arcs_right_slant,ADCs,Heavy_Domains_e,names_list_heavy_e,Light_Domains_e,names_list_light_e,Heavy_Domains_f,names_list_heavy_f,Light_Domains_f,names_list_light_f,Heavy_Domains_g,names_list_heavy_g,Light_Domains_g,names_list_light_g,Heavy_Domains_h,names_list_heavy_h,Light_Domains_h,names_list_light_h,Chem_con)

def render(chains_list,canvas,text_to_image):
    if text_to_image == True:
        canvas.delete("all")
        global all_buttons
        global specificity_colours
        for i in range(len(all_buttons)):
            all_buttons[i].config(fg = "black")
        global canvas_polygons
        canvas_polygons = {}
        global canvas_labels
        canvas_labels= {}
        global Label_lock
        global TYPE_labels
        TYPE_labels = {}
        global NOTE_labels
        NOTE_labels = {}
        global MOD_labels
        MOD_labels   = {}
        global ANTI_labels
        ANTI_labels = {}
        global LENGTH_labels
        LENGTH_labels = {}

        Bonds              = chains_list[0]
        disulphide_bridges = chains_list[1]
        Hinges             = chains_list[2]
        Linkers            = chains_list[3]
        Heavy_Domains_a    = chains_list[4]
        names_Heavy_a      = chains_list[5]
        Light_Domains_a    = chains_list[6]
        names_Light_a      = chains_list[7]
        Heavy_Domains_b    = chains_list[8]
        names_Heavy_b      = chains_list[9]
        Light_Domains_b    = chains_list[10]
        names_Light_b      = chains_list[11]
        Heavy_Domains_c    = chains_list[12]
        names_Heavy_c      = chains_list[13]
        Light_Domains_c    = chains_list[14]
        names_Light_c      = chains_list[15]
        Heavy_Domains_d    = chains_list[16]
        names_Heavy_d      = chains_list[17]
        Light_Domains_d    = chains_list[18]
        names_Light_d      = chains_list[19]
        Heavy_Domains_e    = chains_list[30]
        names_Heavy_e      = chains_list[31]
        Light_Domains_e    = chains_list[32]
        names_Light_e      = chains_list[33]
        Heavy_Domains_f    = chains_list[34]
        names_Heavy_f      = chains_list[35]
        Light_Domains_f    = chains_list[36]
        names_Light_f      = chains_list[37]
        Heavy_Domains_g    = chains_list[38]
        names_Heavy_g      = chains_list[39]
        Light_Domains_g    = chains_list[40]
        names_Light_g      = chains_list[41]
        Heavy_Domains_h    = chains_list[42]
        names_Heavy_h      = chains_list[43]
        Light_Domains_h    = chains_list[44]
        names_Light_h      = chains_list[45]
        Label_Text         = chains_list[20]
        Label_positions    = chains_list[21]
        Domain_Text        = chains_list[22]
        Notes              = chains_list[23]
        Note_positions     = chains_list[24]
        arcs_left          = chains_list[25]
        arcs_right         = chains_list[26]
        arcs_left_slant    = chains_list[27]
        arcs_right_slant   = chains_list[28]
        ADCs               = chains_list[29]
        CCs                = chains_list[46]


    #disulphide_bridge
        if disulphide_bridges != []:
            for i in range(len(disulphide_bridges)):
                domain = canvas.create_line(disulphide_bridges[i], fill='#FF4040', width = 2,tags="disulphide")
                canvas_polygons[domain] = [disulphide_bridges[i], "-disulphide-"]

    #Bonds
        for i in range(len(Bonds)):
            domain = canvas.create_line(Bonds[i], fill=bond_colour, width = 2,tags="bonds")
            canvas_polygons[domain] = [Bonds[i], "-"]
        if arcs_left!=[]:
            for i in range(len(arcs_left)):
                if arcs_left[i][1] == False:
                    domain = canvas.create_arc(arcs_left[i][0], start=90, extent=180, style=tk.ARC, fill=bond_colour, width = 2,tags=("bonds","arcs_left","arcs"))
                    canvas_polygons[domain] = [arcs_left[i][0], "-"]
                elif arcs_left[i][1] == True:
                    domain = canvas.create_arc(arcs_left[i][0], start=90, extent=180, style=tk.ARC, outline=linker_colour, width = 2,tags=("bonds","arcs_left","arcs"))
                    canvas_polygons[domain] = [arcs_left[i][0], "-L-"]
        if arcs_right!=[]:
            for i in range(len(arcs_right)):
                if arcs_right[i][1] == False:
                    domain = canvas.create_arc(arcs_right[i][0], start=270, extent=180, style=tk.ARC, fill=bond_colour, width = 2,tags=("bonds","arcs_right","arcs"))
                    canvas_polygons[domain] = [arcs_right[i][0], "-"]
                elif arcs_right[i][1] == True:
                    domain = canvas.create_arc(arcs_right[i][0], start=270, extent=180, style=tk.ARC, outline=linker_colour, width = 2,tags=("bonds","arcs_right","arcs"))
                    canvas_polygons[domain] = [arcs_right[i][0], "-L-"]
        if arcs_left_slant != []:
            for i in range(len(arcs_left_slant)):
                if arcs_left_slant[i][1] == False:
                    domain = canvas.create_arc(arcs_left_slant[i][0], start=150, extent=120, fill = bond_colour, style=tk.ARC,width=2,tags=("bonds","arcs_left_slant","arcs"))
                    canvas_polygons[domain] = [arcs_left_slant[i][0], "-"]
                elif arcs_left_slant[i][1] == True:
                    domain = canvas.create_arc(arcs_left_slant[i][0], start=150, extent=120, outline = linker_colour, style=tk.ARC,width=2,tags=("bonds","arcs_left_slant","arcs"))
                    canvas_polygons[domain] = [arcs_left_slant[i][0], "-L-"]
        if arcs_right_slant != []:
            for i in range(len(arcs_right_slant)):
                if arcs_right_slant[i][1] == False:
                    domain = canvas.create_arc(arcs_right_slant[i][0], start=270, extent=120, fill = bond_colour, style=tk.ARC,width=2,tags=("bonds","arcs_right_slant","arcs"))
                    canvas_polygons[domain] = [arcs_right_slant[i][0], "-"]
                elif  arcs_right_slant[i][1] == True:
                    domain = canvas.create_arc(arcs_right_slant[i][0], start=270, extent=120, outline = linker_colour, style=tk.ARC,width=2,tags=("bonds","arcs_right_slant","arcs"))
                    canvas_polygons[domain] = [arcs_right_slant[i][0], "-L-"]
        for i in range(len(Linkers)):
            domain = canvas.create_line(Linkers[i], fill=linker_colour, width = 2,tags="bonds")
            canvas_polygons[domain] = [Linkers[i], "-L-"]
        print(Linkers)
        for i in range(len(Hinges)):
            domain = canvas.create_line(Hinges[i], fill=hinge_colour, width = 2,tags=("bonds","hinges"))
            canvas_polygons[domain] = [Hinges[i], "-H-"]

    #A domains
        for i in range(len(Heavy_Domains_a)):
            domain = canvas.create_polygon(Heavy_Domains_a[i], outline='#000000',fill=specificity_colours[0], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_a[i], names_Heavy_a[i]]
        for i in range(len(Light_Domains_a)):
            domain = canvas.create_polygon(Light_Domains_a[i], outline='#000000',fill=specificity_colours[1], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_a[i], names_Light_a[i]]
    #B domains
        for i in range(len(Heavy_Domains_b)):
            domain = canvas.create_polygon(Heavy_Domains_b[i], outline='#000000',fill=specificity_colours[2], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_b[i], names_Heavy_b[i]]
        for i in range(len(Light_Domains_b)):
            domain = canvas.create_polygon(Light_Domains_b[i], outline='#000000',fill=specificity_colours[3], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_b[i], names_Light_b[i]]
    #C domains
        for i in range(len(Heavy_Domains_c)):
            domain = canvas.create_polygon(Heavy_Domains_c[i], outline='#000000',fill=specificity_colours[4], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_c[i], names_Heavy_c[i]]
        for i in range(len(Light_Domains_c)):
            domain = canvas.create_polygon(Light_Domains_c[i], outline='#000000',fill=specificity_colours[5], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_c[i], names_Light_c[i]]
    #D domains
        for i in range(len(Heavy_Domains_d)):
            domain = canvas.create_polygon(Heavy_Domains_d[i], outline='#000000',fill=specificity_colours[6], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_d[i], names_Heavy_d[i]]
        for i in range(len(Light_Domains_d)):
            domain = canvas.create_polygon(Light_Domains_d[i], outline='#000000',fill=specificity_colours[7], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_d[i], names_Light_d[i]]
    #E domains
        for i in range(len(Heavy_Domains_e)):
            domain = canvas.create_polygon(Heavy_Domains_e[i], outline='#000000',fill=specificity_colours[8], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_e[i], names_Heavy_e[i]]
        for i in range(len(Light_Domains_e)):
            domain = canvas.create_polygon(Light_Domains_e[i], outline='#000000',fill=specificity_colours[9], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_e[i], names_Light_e[i]]
    #F domains
        for i in range(len(Heavy_Domains_f)):
            domain = canvas.create_polygon(Heavy_Domains_f[i], outline='#000000',fill=specificity_colours[10], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_f[i], names_Heavy_f[i]]
        for i in range(len(Light_Domains_f)):
            domain = canvas.create_polygon(Light_Domains_f[i], outline='#000000',fill=specificity_colours[11], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_f[i], names_Light_f[i]]
    #G domains
        for i in range(len(Heavy_Domains_g)):
            domain = canvas.create_polygon(Heavy_Domains_g[i], outline='#000000',fill=specificity_colours[12], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_g[i], names_Heavy_g[i]]
        for i in range(len(Light_Domains_g)):
            domain = canvas.create_polygon(Light_Domains_g[i], outline='#000000',fill=specificity_colours[13], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_g[i], names_Light_g[i]]
    #H domains
        for i in range(len(Heavy_Domains_h)):
            domain = canvas.create_polygon(Heavy_Domains_h[i], outline='#000000',fill=specificity_colours[14], width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_h[i], names_Heavy_h[i]]
        for i in range(len(Light_Domains_h)):
            domain = canvas.create_polygon(Light_Domains_h[i], outline='#000000',fill=specificity_colours[15], width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_h[i], names_Light_h[i]]
    #ADCs
        if ADCs != []:
            non_redundant_ADCs = []
            non_redundant_ADCs_sorted = []
            for i in range(len(ADCs)):
                j = sorted(ADCs[i])
                if j not in non_redundant_ADCs_sorted:
                    non_redundant_ADCs.append(ADCs[i])
                    non_redundant_ADCs_sorted.append(j)
            for i in range(len(non_redundant_ADCs)):
                domain = canvas.create_polygon(non_redundant_ADCs[i], outline='#000000',fill=specificity_colours[18], width=2,tags="domain")
                canvas_polygons[domain] = [non_redundant_ADCs[i],  "X"]
    #CCs
        if CCs != []:
            print(CCs)
            non_redundant_CCs = []
            non_redundant_CCs_sorted = []
            for i in range(len(CCs)):
                j = sorted(CCs[i])
                if j not in non_redundant_CCs_sorted:
                    non_redundant_CCs.append(CCs[i])
                    non_redundant_CCs_sorted.append(j)
            for i in range(len(non_redundant_CCs)):
                domain = canvas.create_polygon(non_redundant_CCs[i], outline='#000000',fill=specificity_colours[19], width=2,tags="domain")
                canvas_polygons[domain] = [non_redundant_CCs[i],  "C"]

    #Labels
        if Label_lock == True:
            for i in range(len(Label_positions)):
                x = Label_positions[i][0][0]
                y = Label_positions[i][0][1]
                Domain_Text[i] = re.sub("\_","-", Domain_Text[i])
                Domain_Text[i] = re.sub("nano","",Domain_Text[i])
                label  = canvas.create_text(x,y, text=Domain_Text[i],tags = "label")
                canvas_labels[label] = [[x,y], Domain_Text[i]]


        if Notes != []:
            notes_set = set(Notes)
            setlist = list(notes_set)
            for i in range(len(setlist)):
                if "NOTE:" in setlist[i]:
                    note = canvas.create_text(Note_positions[i],text=setlist[i],tags = "NOTE_labels")
                    NOTE_labels[note] = [Note_positions[i],setlist[i]]
                elif "TYPE:" in setlist[i]:
                    note = canvas.create_text(Note_positions[i],text=setlist[i],tags = "TYPE_labels")
                    TYPE_labels[note] = [Note_positions[i],setlist[i]]
                elif "ANTI:" in setlist[i]:
                    note = canvas.create_text(Note_positions[i],text=setlist[i],tags = "ANTI_labels")
                    ANTI_labels[note] = [Note_positions[i],setlist[i]]
                elif "MOD:" in setlist[i]:
                    note = canvas.create_text(Note_positions[i],text=setlist[i],tags = "MOD_labels")
                    MOD_labels[note] = [Note_positions[i],setlist[i]]
                elif "LENGTH:" in setlist[i]:
                    note = canvas.create_text(Note_positions[i],text=setlist[i],tags = "LENGTH_labels")
                    LENGTH_labels[note] = [Note_positions[i],setlist[i]]





def sequence_render_pipeline(canvas):
    '''
    Gets sequence and automatically tidies it
    '''
    global Bond_lock
    global Delete_lock
    global CustomLabelLock
    global Domain_Primer
    global Domain_Primer_Lock
    global domain_type
    global domain_charge
    global domain_mod
    global all_buttons
    Delete_lock = False
    Bond_lock = ""
    domain_mod = ""
    domain_type = ""
    domain_charge = ""
    Domain_Primer_Lock =""
    CustomLabelLock = ""
    Domain_Primer = []
    for i in range(len(all_buttons)):
        all_buttons[i].config(fg = "black")
    exp = sequence_pipeline(canvas)
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "arrow")
    status_label.config(text="")
    entry=textBox.get("1.0","end-1c")
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates, canvas,True)


############################################
def render_pipeline(canvas):
    '''
    Pipeline takes expression from entry box and generates drawn antibody
    '''
    global Bond_lock
    global Delete_lock
    global CustomLabelLock
    global Domain_Primer
    global Domain_Primer_Lock
    global domain_type
    global domain_charge
    global domain_mod
    global all_buttons
    Delete_lock = False
    Bond_lock = ""
    domain_mod = ""
    domain_type = ""
    domain_charge = ""
    Domain_Primer_Lock =""
    CustomLabelLock = ""
    Domain_Primer = []
    for i in range(len(all_buttons)):
        all_buttons[i].config(fg = "black")
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "arrow")
    status_label.config(text="")
    entry=textBox.get("1.0","end-1c")
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates, canvas,True)

############################################
def prime_domain_button(canvas,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H,domain_name,Light,Heavy):
    global Bond_lock
    global Delete_lock
    global Label_lock
    global Domain_Primer
    global Domain_Primer_Lock
    global all_buttons
    global domain_type
    global domain_mod
    global domain_charge
    global extra_mods
    global Bond_lock
    global Delete_lock
    global Label_lock
    global Domain_Primer_Lock
    global Domain_Primer
    Delete_lock = False
    Bond_lock = ""
    Domain_Primer = []
    global domain_buttons
    global bond_buttons
    for i in domain_buttons:
        i.config(fg = "black")
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    if domain_type != "":
        if Domain_Primer_Lock != domain_name:
            Domain_Primer_Lock = domain_name
            Delete_lock = False
            Bond_lock = ""
            lower_canvas.bind("<Button-1>", mm.place_domain)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
            lower_canvas.config(cursor = "plus")
            status_label.config(text=domain_name)
            Domain_Primer = [righthanded,slant,V,direction,X,mod,interaction,previous_H,domain_name,Light,Heavy]
            if domain_mod != "":
                if  domain_mod == "@" or domain_mod== ">":
                    Domain_Primer[5] == domain_mod

            if domain_charge !="":
                if domain_charge == "+" or domain_charge == "_":
                    domain_charge = re.sub("\_","-", domain_charge)
                    Domain_Primer[8] = re.sub("\+|\-|\_","",Domain_Primer[8])
                    if "V" in Domain_Primer[8]:
                        Domain_Primer[8] = re.sub("\.",domain_charge+".",Domain_Primer[8])
                    else:
                        Domain_Primer[8] += domain_charge
            if extra_mods !="":
                Domain_Primer[8] = re.sub("\!|\*","",Domain_Primer[8])
                if "V" in Domain_Primer[8]:
                    Domain_Primer[8] = re.sub("\.",extra_mods+".",Domain_Primer[8])
                else:
                    Domain_Primer[8] += extra_mods
            if  "nano" in domain_name:
                nanobody_button.config(fg="red")
            elif "VH" in domain_name:
                InsertVHDomainButton.config(fg="red")
            elif "VL" in domain_name:
                InsertVLDomainButton.config(fg="red")
            elif "CH1" in domain_name:
                InsertCH1DomainButton.config(fg="red")
            elif "CH2" in domain_name:
                InsertCH2DomainButton.config(fg="red")
            elif "CH3" in domain_name:
                InsertCH3DomainButton.config(fg="red")
            elif "CH4" in domain_name:
                InsertCH4DomainButton.config(fg="red")
            elif "CL" in domain_name:
                InsertCLDomainButton.config(fg="red")
            elif "X" in domain_name:
                InsertXDomainButton.config(fg="red")
            elif "C" in domain_name:
                InsertCDomainButton.config(fg="red")

        elif Domain_Primer_Lock == domain_name:
            Domain_Primer_Lock = ""
            Bond_lock = ""
            lower_canvas.bind("<Button-1>", mm.select)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.release)
            lower_canvas.config(cursor = "arrow")
            status_label.config(text="")
            Domain_Primer = []
            if domain_type != "":
                Domain_Primer_Lock = ""
                lower_canvas.config(cursor = "arrow")
                lower_canvas.bind("<Button-1>", mm.change_specificity)
                lower_canvas.unbind("<B1-Motion>")
                lower_canvas.unbind("<ButtonRelease-1>")
            #domain_type = ""
            domain_charge = ""
            domain_mod = ""
            extra_mods = ""
    elif domain_type == "":
        domain_type = "a"

        a_button.config(fg="red")
        if Domain_Primer_Lock != domain_name:
            Domain_Primer_Lock = domain_name
            Delete_lock = False
            Bond_lock = ""
            lower_canvas.bind("<Button-1>", mm.place_domain)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
            lower_canvas.config(cursor = "plus")
            status_label.config(text=domain_name)
            Domain_Primer = [righthanded,slant,V,direction,X,mod,interaction,previous_H,domain_name,Light,Heavy]
            if Domain_Primer != []:
                Domain_Primer[8] = re.sub("a|b|c|d|e|f|g|h","",Domain_Primer[8])
                Domain_Primer[8] = re.sub("\.","."+domain_type,Domain_Primer[8])
            if domain_mod != "":
                if  domain_mod == "@" or domain_mod== ">":
                    Domain_Primer[5] == domain_mod

            if domain_charge !="":
                if domain_charge == "+" or domain_charge == "_":
                    domain_charge = re.sub("\_","-", domain_charge)
                    Domain_Primer[8] = re.sub("\+|\-|\_","",Domain_Primer[8])
                    if "V" in Domain_Primer[8]:
                        Domain_Primer[8] = re.sub("\.",domain_charge+".",Domain_Primer[8])
                    else:
                        Domain_Primer[8] += domain_charge
            if extra_mods !="":
                Domain_Primer[8] = re.sub("\!|\*","",Domain_Primer[8])
                if "V" in Domain_Primer[8]:
                    Domain_Primer[8] = re.sub("\.",extra_mods+".",Domain_Primer[8])
                else:
                    Domain_Primer[8] += extra_mods
            if  "nano" in domain_name:
                nanobody_button.config(fg="red")
            elif "VH" in domain_name:
                InsertVHDomainButton.config(fg="red")
            elif "VL" in domain_name:
                InsertVLDomainButton.config(fg="red")
            elif "CH1" in domain_name:
                InsertCH1DomainButton.config(fg="red")
            elif "CH2" in domain_name:
                InsertCH2DomainButton.config(fg="red")
            elif "CH3" in domain_name:
                InsertCH3DomainButton.config(fg="red")
            elif "CH4" in domain_name:
                InsertCH4DomainButton.config(fg="red")
            elif "CL" in domain_name:
                InsertCLDomainButton.config(fg="red")
            elif "X" in domain_name:
                InsertXDomainButton.config(fg="red")
            elif "C" in domain_name:
                InsertCDomainButton.config(fg="red")

        elif Domain_Primer_Lock == domain_name:
            Domain_Primer_Lock = ""
            Bond_lock = ""
            lower_canvas.bind("<Button-1>", mm.select)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.release)
            lower_canvas.config(cursor = "arrow")
            status_label.config(text="")
            Domain_Primer = []
            if domain_type != "":
                Domain_Primer_Lock = ""
                lower_canvas.config(cursor = "arrow")
                lower_canvas.bind("<Button-1>", mm.change_specificity)
                lower_canvas.unbind("<B1-Motion>")
                lower_canvas.unbind("<ButtonRelease-1>")
            #domain_type = ""
            domain_charge = ""
            domain_mod = ""
            extra_mods = ""

############################################
def save_as_png(canvas):
    fileName = "AbYdraw_export"
    # save postscipt image

    #canvas.postscript(file = fileName + '.eps')
    # use PIL to convert to PNG


    #f = filedialog.asksaveasfile(mode='wb', defaultextension=".png")
    #if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
    #    return
    #f.save(img)
    #f.close()
    #
    file = filedialog.asksaveasfile(mode='w', defaultextension=".eps", filetypes=(("EPS file", "*.eps"),("All Files", "*.*") ))
    if file:
        abs_path = os.path.abspath(file.name)
        eps = canvas.postscript(file = abs_path,colormode='color',height=1000)
        file.close()
        #eps.save(abs_path, 'eps', lossless = True) # saves the image to the input file name.
        #s.remove(str(fileName + '.eps'))

###########################################
def get_min_max_coordinates(domain_coordinates):
    xs= []
    ys= []
    for f in range(len(domain_coordinates)):
        if f%2==0:
            xs.append(domain_coordinates[f])
        elif f %2 !=0:
            ys.append(domain_coordinates[f])
    x1 = min(xs)
    x2 = max(xs)
    y1 = min(ys)
    y2 = max(ys)
    avx= (x1+x2)/2
    avy= (y1+y2)/2
    return(x1,x2,y1,y2,avx,avy)
#########################################
def Get_Template_File(canvas):
    '''
    Take in entry and export template file
    '''
    global Bond_lock
    global Delete_lock
    global textBox
    Delete_lock = False
    Bond_lock = ""
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "arrow")
    status_label.config(text="")
    entry=textBox.get("1.0","end-1c")
    Template_File = []
    Sequences = []
    chains_list = Get_dictionaries(entry)
    VHa_chain       = chains_list[0]
    VLa_chain       = chains_list[1]
    VHb_chain       = chains_list[2]
    VLb_chain       = chains_list[3]
    chain_count     = chains_list[4]
    fragment1       = chains_list[5]
    fragment2       = chains_list[6]
    fragment3       = chains_list[7]
    fragment4       = chains_list[8]
    VHa_chain_list  = list(VHa_chain.keys())
    VLa_chain_list  = list(VLa_chain.keys())
    VHb_chain_list  = list(VHb_chain.keys())
    VLb_chain_list  = list(VLb_chain.keys())
    fragment1_list  = list(fragment1.keys())
    fragment2_list  = list(fragment2.keys())
    fragment3_list  = list(fragment3.keys())
    fragment4_list  = list(fragment4.keys())


    def Get_Domain_stats_printout(string,dict,indexing):
        #print(string)
        #print(dict)
        Type = []
        Chain_Class_Strings = []
        NGlycos_Strings = []
        Heavy_CysPositions = []
        Heavy_Disulphides_Intra = []
        Light_CysPositions = []
        Light_Disulphides_Intra = []
        DisulphidesInter = []
        HVJC_Germline_Strings = []
        H_chain_Ranges_Strings = []
        H_Mutations = []
        LVJC_Germline_Strings = []
        L_chain_Ranges_Strings = []
        L_Mutations = []
        H_CDR_Strings = []
        L_CDR_Strings =[]
        all = [Type,Chain_Class_Strings,NGlycos_Strings,Heavy_CysPositions,Heavy_Disulphides_Intra,Light_CysPositions,Light_Disulphides_Intra,DisulphidesInter,HVJC_Germline_Strings,H_chain_Ranges_Strings,H_Mutations,LVJC_Germline_Strings,L_chain_Ranges_Strings,L_Mutations,H_CDR_Strings,L_CDR_Strings, [], []]
        index = str("["+str(indexing)+"]")
        #dict_keys_list = list(dict.keys())
        #disulphides_connected = []
        #disulphides_inter = 0
        #disulphides_intra = 0
        #for i in range(len(dict_keys_list)):
        #    if dict.get(dict_keys_list[i])[1] > 0:
        #        disulphides_connected.append(dict.get(dict_keys_list[i])[0])[1]
        #for i in range(len(disulphides_connected)):
        #    found= False
        #    for j in range(len(dict_keys_list)):
        #        if disulphides_connected == dict.get(dict_keys_list[j])[0]:
        #            disulphides_intra += 0.5
        #            found = True
        #    if found == False:
        #        disulphides_inter +=1
        #if disulphides_intra > 0:
        #    Heavy_Disulphides_Intra.append(str("HeavyCysPositions"+index+":"))
        #    Heavy_Disulphides_Intra.append(str("HeavyPotentialNGlycos"+index+":"))
        #    if "L" in str(string):
        #        Light_Disulphides_Intra.append(str("LightCysPositions"+index+":"))
        #        Light_Disulphides_Intra.append(str("LightPotentialNGlycos"+index+":"))

        #if disuphides_inter > 0:
        #    DisulphidesInter.append(str("DisulfidesInter"+index+":"))

        if "X" in str(string):
            Type.append(str("Type"+index+": OTHER"))
        else:
            Type.append(str("Type"+index+":"))
        if "VH" in str(string):
            #Chain_class
            NGlycos_Strings.append(str("HeavyPotentialNGlycos"+index+":"))
            NGlycos_Strings.append(str("HeavyConfirmedNGlycos"+index+":"))
            #Cys_pos_and_disulphides

            #Germline
            HVJC_Germline_Strings.append(str("HVGermline"+index+":"))
            HVJC_Germline_Strings.append(str("HJGermline"+index+":"))
            #Ranges
            H_chain_Ranges_Strings.append(str("VHRange"+index+":"))
            #CDRs
            H_CDR_Strings.append(str("CDRKabatH1"+index+":"))
            H_CDR_Strings.append(str("CDRKabatH2"+index+":"))
            H_CDR_Strings.append(str("CDRKabatH3"+index+":"))
        if "VL" in str(string):
            #Chain_class
            NGlycos_Strings.append(str("LightPotentialNGlycos"+index+":"))
            NGlycos_Strings.append(str("LightConfirmedNGlycos"+index+":"))
            #Cys_pos_and_disulphides
            #Germline
            LVJC_Germline_Strings.append(str("LVGermline"+index+":"))
            LVJC_Germline_Strings.append(str("LJGermline"+index+":"))
            #Ranges
            L_chain_Ranges_Strings.append(str("VLRange"+index+":"))
            #CDRs
            L_CDR_Strings.append(str("CDRKabatL1"+index+":"))
            L_CDR_Strings.append(str("CDRKabatL2"+index+":"))
            L_CDR_Strings.append(str("CDRKabatL3"+index+":"))
        if "CH1" in str(string):
            HVJC_Germline_Strings.append(str("HCGermline"+index+":"))
            H_chain_Ranges_Strings.append(str("CH1Range"+index+":"))
        if "H" in str(string):
            H_chain_Ranges_Strings.append(str("HingeRange"+index+":"))
        if "CH2" in str(string):
            H_chain_Ranges_Strings.append(str("CH2Range"+index+":"))
        if "CH3" in str(string):
            H_chain_Ranges_Strings.append(str("CH3Range"+index+":"))
        if "CH4" in str(string):
            H_chain_Ranges_Strings.append(str("CH4Range"+index+":"))
        if "CH1" in str(string) or "CH2" in str(string) or "CH3" in str(string) or "CH4" in str(string):
            H_chain_Ranges_Strings.append(str("CHSRange"+index+":"))
        if "CL" in str(string):
            LVJC_Germline_Strings.append(str("LCGermline"+index+":"))
            L_chain_Ranges_Strings.append(str("CLRange"+index+":"))
        if "H@" in str(string):
            H_Mutations.append(str("MutationH"+index+":"+"          heterodimer formation hole "))
        if "H>" in str(string):
            L_Mutations.append(str("MutationH"+index+":"+"          heterodimer formation knob "))
        if "L@" in str(string):
            L_Mutations.append(str("MutationH"+index+":"+"          heterodimer formation hole "))
        if "L>" in str(string):
            H_Mutations.append(str("MutationH"+index+":"+"          heterodimer formation knob "))
        Template_File.append("")
        Template_File.append("")
        for j in range(len(all)):
            if all[j] != []:
                for k in range(len(all[j])):
                    Template_File.append(str(all[j][k]))

    def Get_Linker_stats_printout(j):
        Template_File.append("")
        Template_File.append("")
        Template_File.append(str("Linker["+str(j)+","+str(j+1)+"]:"))
        Template_File.append("")
        Template_File.append("")

    def save_txt_file(Template_File):
        to_save = ""
        for i in Template_File:
            to_save += (str(i)+"\n")
        f = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
        if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
            return
        f.write(to_save)
        f.close()



    if (VHa_chain_list  == ['VHa[0]', 'CH1[1]', 'H[2]', 'CH2[3]', 'CH3[4]'] and VLa_chain_list == ['VLa[0]', 'CL[1]']) or (VHa_chain_list  == ['VHa[0]', 'CH1[1]', 'H[2]'] and VLa_chain_list == ['VLa[0]', 'CL[1]']) or (VHa_chain_list == ['VHb[0]', 'CH1[1]', 'H[2]', 'CH2[3]', 'CH3[4]'] and VLb_chain_list == ['VLb[0]', 'CL[1]']) or (VHa_chain_list == ['VHb[0]', 'CH1[1]', 'H[2]'] and VLa_chain_list == ['VLb[0]', 'CL[1]']):
        for i in range(len(VLa_chain_list)):
            VHa_chain[VLa_chain_list[i]] = VLa_chain.get(VLa_chain_list[i])
            VHa_chain_list.append(VLa_chain_list[i])
        VLa_chain = {}
        VLa_chain_list = []
    if (VHb_chain_list  == ['VHb[0]', 'CH1[1]', 'H[2]', 'CH2[3]', 'CH3[4]'] and VLb_chain_list == ['VLb[0]', 'CL[1]']) or (VHb_chain_list  == ['VHb[0]', 'CH1[1]', 'H[2]'] and VLb_chain_list == ['VLb[0]', 'CL[1]']) or (VHb_chain_list  == ['VHa[0]', 'CH1[1]', 'H[2]', 'CH2[3]', 'CH3[4]'] and VLb_chain_list == ['VLa[0]', 'CL[1]']) or (VHb_chain_list  == ['VHa[0]', 'CH1[1]', 'H[2]'] and VLb_chain_list == ['VLa[0]', 'CL[1]']):
        for i in range(len(VLb_chain_list)):
            VHb_chain[VLb_chain_list[i]] = VLb_chain.get(VLb_chain_list[i])
            VHb_chain_list.append(VLb_chain_list[i])
        VLb_chain = {}
        VLb_chain_list = []

    all_chains_dict = [VHa_chain,VLa_chain,VHb_chain,VLb_chain,fragment1,fragment2,fragment3,fragment4]
    all_chains_keyslist= [VHa_chain_list,VLa_chain_list,VHb_chain_list,VLb_chain_list,fragment1_list,fragment2_list,fragment3_list,fragment4_list]
    all_chains_keyslist = [x for x in all_chains_keyslist if x != []]
    all_chains_dict = [x for x in all_chains_dict if x != {}]
    fusions        = []
    fusions_dicts  = []
    fusion_indexes = []
    fusion_counter = 1
    for i in range(len(all_chains_keyslist)):
        chain_fusions_indexes     = []
        chain_fusions            = []
        chain_fusions_dict        =[]
        current_fusion_index     = []
        current_fusion     = []
        current_fusion_dict={}
        for j in range(len(all_chains_keyslist[i])):
            if j+1 != len(all_chains_keyslist[i]):
                if "Linker[" not in str(all_chains_keyslist[i][j]):
                    current_fusion.append(all_chains_keyslist[i][j])
                    current_fusion_dict[all_chains_keyslist[i][j]] = all_chains_dict[i].get(all_chains_keyslist[i][j])
                elif "Linker[" in str(all_chains_keyslist[i][j]):
                    chain_fusions.append(current_fusion)
                    chain_fusions_dict.append(current_fusion_dict)
                    current_fusion_index.append(fusion_counter)
                    fusion_counter +=1
                    current_fusion = []
                    current_fusion_dict = {}
            elif j+1 == len(all_chains_keyslist[i]):
                current_fusion.append(all_chains_keyslist[i][j])
                current_fusion_dict[all_chains_keyslist[i][j]] = all_chains_dict[i].get(all_chains_keyslist[i][j])
                chain_fusions.append(current_fusion)
                chain_fusions_dict.append(current_fusion_dict)
                current_fusion_index.append(fusion_counter)
                chain_fusions_indexes.append(current_fusion_index)
                fusion_counter +=1
                current_fusion = []
                current_fusion_dict = {}
                current_fusion_index = []
        fusions.append(chain_fusions)
        fusions_dicts.append(chain_fusions_dict)
        fusion_indexes.append(chain_fusions_indexes)


    fusions_dicts = [x for x in fusions_dicts if x != []]
    fusion_indexes= [x for x in fusion_indexes if x != []]

    paired_V_domains = []
    already_done = []
    V_domains = {}
    counter = 1
    for i in range(len(fusions)):
        for j in range(len(fusions[i])):
            mini_counter = 1
            for k in range(len(fusions[i][j])):
                if "V" in str(fusions[i][j][k]):
                    if mini_counter == 1:
                        V_domains[counter] = fusions_dicts[i][j].get(fusions[i][j][k])[0]
                        mini_counter +=1
                    elif mini_counter >1:
                        paired_V_domains.append(counter)
                        del V_domains[counter]
            counter+=1
    V_domains_keys = list(V_domains.keys())
    for i in range(len(V_domains_keys)):
        interactor = V_domains.get(V_domains_keys[i])[0]
        for j in range(len(V_domains_keys)):
            current_interactor = V_domains.get(V_domains_keys[j])[1]
            if interactor == current_interactor:
                to_append = str("Antigen["+str(i+1)+","+str(j+1)+"]:")
                to_done = str("Antigen["+str(j+1)+","+str(i+1)+"]:")
                if to_append not in already_done:
                    paired_V_domains.append(to_append)
                    already_done.append(to_done)
    for i in paired_V_domains:
        Template_File.append(i)


    counter = 1
    for i in range(len(fusions_dicts)):
        for j in range(len(fusions_dicts[i])):
            listing = list(fusions_dicts[i][j].keys())
            listier = re.sub("\[[0-9]\]|\[|\]|a|b|c|d|e|f|g|h|'|,","",str(listing))
            index = str(counter)
            counter +=1
            Template_File.append(str("Domains["+index+"]   "+listier))


    for i in range(len(fusion_indexes)):
        if fusion_indexes[i] != []:
            listier = re.sub("\[|\]|'|,","",str(fusion_indexes[i]))
            Template_File.append(str("Fusion    "+str(listier) ))




    counter = 1
    for i in range(len(fusions_dicts)):
        sequences_to_append = []
        for j in range(len(fusions_dicts[i])):
            listing = list(fusions_dicts[i][j].keys())
            listier = re.sub("\[[0-9]\]|\[|\]|'|,","",str(listing))
            index = str(counter)
            sequences_to_append.append(counter)
            Get_Domain_stats_printout(listing,fusions_dicts,index)
            if j+1 != len(fusions_dicts[i]):
                Get_Linker_stats_printout(counter)
            counter +=1
        print(str(listing))
        if "VH" in str(listing) and "VL" in str(listing):
            Sequences.append("")
            Sequences.append(str("Chain"+str(sequences_to_append)+"Heavy:"))
            Sequences.append("")
            Sequences.append(str("Chain"+str(sequences_to_append)+"Light:"))
            Sequences.append("")
        else:
            Sequences.append("")
            Sequences.append(str("Chain"+str(sequences_to_append)+":"))
            Sequences.append("")



    for i in Sequences:
        Template_File.append(i)
    String = textBox.get("1.0","end")
    Template_File.append(str("AbML String: "+String))

    save_txt_file(Template_File)




############################################
def sequence_pipeline(canvas):
    '''
    Take drawing from cavnas and generate expression to be displayed on entry box
    '''
    global all_buttons
    for i in range(len(all_buttons)):
        all_buttons[i].config(fg="black")
    lower_canvas.config(cursor = "arrow")
    Bond_lock = ""
    Delete_lock = False
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    polygons_keyslist = list(canvas_polygons.keys())
    TYPE_keyslist = list(TYPE_labels.keys())
    NOTE_keyslist = list(NOTE_labels.keys())
    MOD_keyslist  = list(MOD_labels.keys())
    ANTI_keyslist = list(ANTI_labels.keys())
    LENGTH_keyslist=list(LENGTH_labels.keys())
    domains_list = canvas.find_withtag("domain")
    domains_dict = {}
    for i in range(len(domains_list)):
        for j in range(len(polygons_keyslist)):
            if domains_list[i] == polygons_keyslist[j]:
                domains_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    bonds_list = canvas.find_withtag("bonds")
    bonds_dict = {}
    for i in range(len(bonds_list)):
        for j in range(len(polygons_keyslist)):
            if bonds_list[i] == polygons_keyslist[j]:
                bonds_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    disulphides_list = canvas.find_withtag("disulphide")
    disulphides_dict = {}
    for i in range(len(disulphides_list)):
        for j in range(len(polygons_keyslist)):
            if disulphides_list[i] == polygons_keyslist[j]:
                disulphides_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    type_list = canvas.find_withtag("TYPE_labels")
    type_dict = {}
    for i in range(len(type_list)):
        for j in range(len(TYPE_keyslist)):
            if type_list[i] == TYPE_keyslist[j]:
                type_dict[j] = TYPE_labels.get(TYPE_keyslist[j])
    note_list = canvas.find_withtag("NOTE_labels")
    note_dict = {}
    for i in range(len(note_list)):
        for j in range(len(NOTE_keyslist)):
            if note_list[i] == NOTE_keyslist[j]:
                note_dict[j] = NOTE_labels.get(NOTE_keyslist[j])
    mod_list = canvas.find_withtag("MOD_labels")
    mod_dict = {}
    for i in range(len(mod_list)):
        for j in range(len(MOD_keyslist)):
            if mod_list[i] == MOD_keyslist[j]:
                mod_dict[j] = MOD_labels.get(MOD_keyslist[j])
    anti_list = canvas.find_withtag("ANTI_labels")
    anti_dict = {}
    for i in range(len(anti_list)):
        for j in range(len(ANTI_keyslist)):
            if anti_list[i] == ANTI_keyslist[j]:
                anti_dict[j] = ANTI_labels.get(ANTI_keyslist[j])
    length_list = canvas.find_withtag("LENGTH_labels")
    length_dict = {}
    for i in range(len(length_list)):
        for j in range(len(LENGTH_keyslist)):
            if length_list[i] == LENGTH_keyslist[j]:
                length_dict[j] = LENGTH_labels.get(LENGTH_keyslist[j])
    arcs_list = canvas.find_withtag("arcs")
    arcs_dict = {}
    for i in range(len(arcs_list)):
        for j in range(len(polygons_keyslist)):
            if arcs_list[i] == polygons_keyslist[j]:
                arcs_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    hinges_list = canvas.find_withtag("hinges")
    hinges_dict = {}
    for i in range(len(hinges_list)):
        for j in range(len(polygons_keyslist)):
            if hinges_list[i] == polygons_keyslist[j]:
                hinges_dict[j] = canvas_polygons.get(polygons_keyslist[j])

    def arc_checker(key):
        bondx1 = bonds_dict.get(key)[0][0]
        bondy1 = bonds_dict.get(key)[0][1]
        bondx2 = bonds_dict.get(key)[0][2]
        bondy2 = bonds_dict.get(key)[0][3]

        if bondx2-bondx1 == 100: #arcs_left
            bondx1 = bondx1+50
            bondx2 = bondx2-50
            bondy1 = bondy1+20
        elif bondx1-bondx2 == 100:
            bondx1 = bondx1-50
            bondx2 = bondx2+50
            bondy1 = bondy1+20
        elif bondx2>bondx1:
            bondx2 = bondx2-50
            bondx1,bondx2 = bondx2,bondx1
        elif bondx1>bondx2:
            bondx2 = bondx2+50
            bondx1,bondx2 = bondx2,bondx1

        upside_down = False
        if bondy1<bondy2:
            for x in range(len(domains_keyslist)):
                min_max = get_min_max_coordinates(domains_dict.get(domains_keyslist[i])[0])
                domainx1 = min_max[0]
                domainx2 = min_max[1]
                domainy1 = min_max[2]
                domainy2 = min_max[3]
                domainy = min_max[5]
                if  (domainx1<= bondx1 <= domainx2) and (domainy<=bondy1<domainy2):
                    upside_down = True
        if upside_down == False:
            bondy1,bondy2=bondy2,bondy1


        bondy2= bondy2+5

        return(bondx1,bondy1,bondx2,bondy2)
    print("DOMAINs", domains_dict)
    print("BONDS", bonds_dict)
    print("Disulphides", disulphides_dict)
    print("Types", type_dict)
    print("Notes", note_dict)
    print("anti", anti_dict)
    print("mod", mod_dict)
    print("length", length_dict)
    print("ARCS", arcs_dict)
    chains=[]
    current_chain_str = []
    current_chain_coords_lists = []
    domains_keyslist = list(domains_dict.keys())
    bonds_keyslist = list(bonds_dict.keys())
    disulphides_keyslist = list(disulphides_dict.keys())
    type_keyslist = list(type_dict.keys())
    note_keyslist = list(note_dict.keys())
    mod_keyslist  = list(mod_dict.keys())
    anti_keyslist = list(anti_dict.keys())
    length_keyslist= list(length_dict.keys())
    arcs_keyslist = list(arcs_dict.keys())
    hinges_keyslist= list(hinges_dict.keys())
    for i in range(len(domains_keyslist)):
        start_found = True
        continuation_found = False
        min_max = get_min_max_coordinates(domains_dict.get(domains_keyslist[i])[0])
        domainx1 = min_max[0]
        domainx2 = min_max[1]
        domainy1 = min_max[2]
        domainy2 = min_max[3]
        for j in range(len(bonds_keyslist)):
            bondx2 = bonds_dict.get(bonds_keyslist[j])[0][2]
            bondy2 = bonds_dict.get(bonds_keyslist[j])[0][3]
            if arcs_keyslist != []:
                if bonds_keyslist[j] in arcs_keyslist:
                    bondx2 = arc_checker(bonds_keyslist[j])[2]
                    bondy2 = arc_checker(bonds_keyslist[j])[3]
            if domainx1 < bondx2 < domainx2 and domainy1 < bondy2 < domainy2:
                start_found = False

        if start_found == True:
            chains.append(domains_keyslist[i])
            start_found = False

    strings         = []
    full_chains     = []
    #full_directions = []



##get raw chains in format
    for i in range(len(chains)):
        start = chains[i]
        string = []
        full_chain = []
        continuation_found = True
        while continuation_found == True:
            if full_chain == []:
                connection_found = False
                full_chain.append(chains[i])
                string.append(domains_dict.get(chains[i])[1])
                #directions.append(domains_dict.get(chains[i])[2])
                min_max = get_min_max_coordinates(domains_dict.get(full_chain[-1])[0])
                domainx1 = min_max[0]
                domainx2 = min_max[1]
                domainy1 = min_max[2]
                domainy2 = min_max[3]
                for j in range(len(bonds_keyslist)):
                    bondx1 = bonds_dict.get(bonds_keyslist[j])[0][0]
                    bondy1 = bonds_dict.get(bonds_keyslist[j])[0][1]
                    if arcs_keyslist != []:
                        if bonds_keyslist[j] in arcs_keyslist:
                            bondx1 = arc_checker(bonds_keyslist[j])[0]
                            bondy1 = arc_checker(bonds_keyslist[j])[1]
                    if domainx1 < bondx1 < domainx2 and domainy1 < bondy1 < domainy2:
                        connection_found = True
                        full_chain.append(bonds_keyslist[j])
                        string.append(bonds_dict.get(bonds_keyslist[j])[1])
                        break
                if connection_found == True:
                    continuation_found = True
                else:
                    continuation_found = False

            elif full_chain !=[]:
                connection_found = False
                ##find continuing bond
                bondx2 = bonds_dict.get(full_chain[-1])[0][2]
                bondy2 = bonds_dict.get(full_chain[-1])[0][3]
                if arcs_keyslist != []:
                    if full_chain[-1] in arcs_keyslist:
                        bondx2 = arc_checker(full_chain[-1])[2]
                        bondy2 = arc_checker(full_chain[-1])[3]


                for j in range(len(domains_dict)):
                    min_max = get_min_max_coordinates(domains_dict.get(domains_keyslist[j])[0])
                    domainx1 = min_max[0]
                    domainx2 = min_max[1]
                    domainy1 = min_max[2]
                    domainy2 = min_max[3]
                    if domainx1 < bondx2 < domainx2 and domainy1 < bondy2 < domainy2:
                        continuation_found = True
                        full_chain.append(domains_keyslist[j])
                        string.append(domains_dict.get(domains_keyslist[j])[1])
                        ##Check for another bond
                        min_max = get_min_max_coordinates(domains_dict.get(full_chain[-1])[0])
                        domainx1 = min_max[0]
                        domainx2 = min_max[1]
                        domainy1 = min_max[2]
                        domainy2 = min_max[3]
                        for n in range(len(bonds_keyslist)):
                            bondx1 = bonds_dict.get(bonds_keyslist[n])[0][0]
                            bondy1 = bonds_dict.get(bonds_keyslist[n])[0][1]
                            if arcs_keyslist != []:
                                if bonds_keyslist[n] in arcs_keyslist:
                                    bondx1 = arc_checker(bonds_keyslist[n])[0]
                                    bondy1 = arc_checker(bonds_keyslist[n])[1]
                            if domainx1 < bondx1 < domainx2 and domainy1 < bondy1 < domainy2:
                                connection_found = True
                                full_chain.append(bonds_keyslist[n])
                                string.append(bonds_dict.get(bonds_keyslist[n])[1])
                                break


                if connection_found == True:
                    continuation_found = True
                else:
                    continuation_found = False
                    for j in range(len(hinges_keyslist)): #check for extra bonds
                        try:
                            bondx2 = bonds_dict.get(full_chain[-1])[0][2]
                            bondy2 = bonds_dict.get(full_chain[-1])[0][3]
                            if arcs_keyslist != []:
                                if full_chain[-1] in arcs_keyslist:
                                    bondx2 = arc_checker(full_chain[-1])[2]
                                    bondy2 = arc_checker(full_chain[-1])[3]
                            hingex1 = hinges_dict.get(hinges_keyslist[j])[0][0]
                            hingey1 = hinges_dict.get(hinges_keyslist[j])[0][1]
                            print(bondx2,hingex1,bondx2-hingex1)
                            print(bondy2,hingey1,bondy2-hingey1 )
                            if -25<=(bondx2-hingex1)<=25 and -25<=(bondy2-hingey1)<=25:
                                full_chain.append(hinges_keyslist[j])
                                print(hinges_dict.get(hinges_keyslist[j]))
                                string.append(hinges_dict.get(hinges_keyslist[j])[1])
                        except TypeError:
                            pass




        #full_directions.append(directions)
        full_chains.append(full_chain)
        strings.append(string)

##number chains
    counter = 1
    assigned_numbers = {}
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            assigned_keyslist = list(assigned_numbers.keys())
            if str(strings[i][j]) != "-" and  str(strings[i][j]) != "-L-":
                if "-H-" not in strings[i][j]:
                    index = full_chains[i][j]
                    coordinates = domains_dict.get(index)[0]

                    assigned_match = False
                    assigned_keyslist_check = ""
                    for k in range(len(assigned_keyslist)):
                        assigned_coordinates = assigned_numbers.get(assigned_keyslist[k])
                        #print("OK THIS IS A RAID", coordinates, assigned_coordinates)
                        if coordinates == assigned_coordinates:
                            assigned_match = True
                            assigned_keyslist_check = assigned_keyslist[k]
                            print(assigned_keyslist[k],coordinates, assigned_coordinates)
                    if assigned_match == True:
                        if "X" in strings[i][j]:
                            strings[i][j] = str("-X("+str(assigned_keyslist_check)+")-")
                        elif strings[i][j] == "C":
                            strings[i][j] = str("-C("+str(assigned_keyslist_check)+")-")
                    elif assigned_match == False:
                        if "-" not in strings[i][j]:
                            strings[i][j] += str("("+str(counter)+")")
                        elif "-X-" in strings[i][j]:
                            strings[i][j] = str("-X("+str(counter)+")-")
                        elif  strings[i][j] == "C":
                            strings[i][j] = str("-C("+str(counter)+")-")
                        assigned_numbers[counter] = coordinates
                        counter += 1
                elif "-H-" in strings[i][j]:
                    strings[i][j] = str("-H("+str(counter)+")-")
                    counter += 1


    print(strings)
##Pair chains based on closeness
    paired = []
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            print(i,j)
            if ":" not in str(strings[i][j]) and "-" not in str(strings[i][j]) and "nano" not in str(strings[i][j]):
                number =  re.findall("\((.*?)\)", str(strings[i][j]))
                number =  int(re.sub("\[|\'|\]","", str(number)))

                if number not in paired and ("X[" not in str(strings[i][j]) and str(strings[i][j]) != "C"):
                    domain_name = re.sub("\((.*?)\)","",str(strings[i][j]))
                    index = full_chains[i][j]
                    min_max = get_min_max_coordinates((domains_dict.get(index)[0]))
                    d1x1 = min_max[0]
                    d1x2 = min_max[1]
                    d1y1 = min_max[2]
                    d1y2 = min_max[3]
                    if "V" in str(domains_dict.get(index)[1]):
                        testx1 = (domains_dict.get(index)[0][4])
                        testx2 = (domains_dict.get(index)[0][8])
                        if testx1 > testx2:
                            d1x1 = min_max[0]-30
                        else:
                            d1x2 = min_max[1]+30
                    else:
                        d1x1 = min_max[0]-30
                        d1x2 = min_max[1]+30
                        d1y1 = min_max[2]-5
                        d1y2 = min_max[3]+5
                    #print(number)
                    #print(d1x1,d1x2,d1y1,d1y2)
                    #Search for overlapping matching domains
                    for f in range(len(domains_keyslist)):
                        if domains_keyslist[f] != index:
                            for a in range(len(full_chains)):
                                for b in range(len(full_chains[a])):
                                    if int(full_chains[a][b]) == int(domains_keyslist[f]) and ":" not in strings[a][b]:
                                        paired_domain = strings[a][b]
                                        paired_number =  re.findall("\((.*?)\)", str(paired_domain))
                                        paired_number =  str(re.sub("\[|\'|\]","", str(paired_number)))
                                        paired_name   = re.sub("\((.*?)\)", "", str(paired_domain))
                                        in_paired = False
                                        for x in range(len(paired)):
                                            if paired[x] == int(number):
                                                in_paired = True
                                        if in_paired == False:
                                            min_max = get_min_max_coordinates(domains_dict.get(domains_keyslist[f])[0])
                                            d2x1 = min_max[0]
                                            d2x2 = min_max[1]
                                            d2y1 = min_max[2]
                                            d2y2 = min_max[3]
                                            if "V" in str(domains_dict.get(index)[1]):
                                                testx1 = (domains_dict.get(index)[0][4])
                                                testx2 = (domains_dict.get(index)[0][8])
                                                if testx1 > testx2:
                                                    d2x1 = min_max[0]-30
                                                else:
                                                    d2x2 = min_max[1]+30
                                            else:
                                                d2x1 = min_max[0]-30
                                                d2x2 = min_max[1]+30
                                                d2y1 = min_max[2]-5
                                                d2y2 = min_max[3]+5
                                            combinations_to_try = [[d2x1,((d2y1+d2y2)/2)],[d2x2,((d2y1+d2y2)/2)],[d2x1,d2y1],[d2x2,d2y1],[d2x2,d2y2],[d2x1,d2y2]]
                                            for g in combinations_to_try:
                                                if d1x1 < g[0] < d1x2 and d1y1 < g[1] < d1y2:
                                                    if ("VH" in str(strings[i][j]) and "VL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("VL" in str(strings[i][j]) and "VH" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CL" in str(strings[i][j]) and "CH1" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH1" in str(strings[i][j]) and "CL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH2" in str(strings[i][j]) and "CH2" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH3" in str(strings[i][j]) and "CH3" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH4" in str(strings[i][j]) and "CH4" in str(domains_dict.get(domains_keyslist[f])[1])) or ("-H-" == str(strings[i][j]) and "-H-" == str(domains_dict.get(domains_keyslist[f])[1])) or ("X" in str(strings[i][j]) and "X" in str(domains_dict.get(domains_keyslist[f])[1])) or ( str(strings[i][j]) == "C" and  str(domains_dict.get(domains_keyslist[f])[1]) == "C"):
                                                        disulphide_count = 0
                                                        for y in range(len(disulphides_keyslist)):
                                                            #print("Looking for those disulphides")
                                                            disulphx1 = disulphides_dict.get(disulphides_keyslist[y])[0][0]
                                                            disulphy1 = disulphides_dict.get(disulphides_keyslist[y])[0][1]
                                                            disulphx2 = disulphides_dict.get(disulphides_keyslist[y])[0][2]
                                                            disulphy2 = disulphides_dict.get(disulphides_keyslist[y])[0][3]
                                                            #print("disulphx1",d1x1, disulphx1, d1x2 ,"disulphy1", d1y1, disulphy1 , d1y2,"disulphx2" , d2x1 , disulphx2 , d2x2 ,"disulphy2", d2y1, disulphy2 , d2y2)
                                                            #print("disulphx2",d1x1, disulphx2, d1x2 ,"disulphy2", d1y1, disulphy2 , d1y2,"disulphx1" , d2x1 , disulphx1 , d2x2 ,"disulphy1", d2y1, disulphy1 , d2y2)

                                                            if ((d1x1 <= disulphx2 <= d1x2 and d1y1 <= disulphy2 <= d1y2) and (d2x1 <= disulphx1 <= d2x2 and d2y1 <= disulphy1 <= d2y2)) or ((d1x1 <= disulphx1 <= d1x2 and d1y1 <= disulphy1 <= d1y2) and (d2x1 <= disulphx2 <= d2x2 and d2y1 <= disulphy2 <= d2y2)):
                                                                disulphide_count += 1

                                                        if disulphide_count == 0:
                                                            strings[i][j] = str(domain_name+"("+str(number)+":"+str(paired_number)+")")
                                                            strings[a][b] = str(paired_name+"("+str(paired_number)+":"+str(number)+")")
                                                        elif disulphide_count > 0:
                                                            strings[i][j] = str(domain_name+"("+str(number)+":"+str(paired_number)+")"+"{"+str(disulphide_count)+"}")
                                                            strings[a][b] = str(paired_name+"("+str(paired_number)+":"+str(number)+")"+"{"+str(disulphide_count)+"}")
                                                        paired.append(int(number))
                                                        paired.append(int(paired_number))

                elif ("X" in str(strings[i][j]) or str(strings[i][j])) == "C" :
                    domain_name = re.sub("\((.*?)\)","",str(strings[i][j]))
                    paired_X_domains = []
                    index = full_chains[i][j]
                    min_max = get_min_max_coordinates((domains_dict.get(index)[0]))
                    d1x1 = min_max[0]-30
                    d1x2 = min_max[1]+30
                    d1y1 = min_max[2]-20
                    d1y2 = min_max[3]+20
                    #Search for overlapping matching domains
                    for f in range(len(domains_keyslist)):
                        if domains_keyslist[f] != index:
                            for a in range(len(full_chains)):
                                for b in range(len(full_chains[a])):
                                    if int(full_chains[a][b]) == int(domains_keyslist[f]) and ":" not in strings[a][b]:
                                        paired_domain = strings[a][b]
                                        paired_number =  re.findall("\((.*?)\)", str(paired_domain))
                                        paired_number =  str(re.sub("\[|\'|\]","", str(paired_number)))
                                        paired_name   = re.sub("\((.*?)\)", "", str(paired_domain))

                                        min_max = get_min_max_coordinates(domains_dict.get(domains_keyslist[f])[0])
                                        d2x1 = min_max[0]-30
                                        d2x2 = min_max[1]+30
                                        d2y1 = min_max[2]-20
                                        d2y2 = min_max[3]+20
                                        combinations_to_try = [[d2x1,((d2y1+d2y2)/2)],[d2x2,((d2y1+d2y2)/2)],[d2x1,d2y1],[d2x2,d2y1],[d2x2,d2y2],[d2x1,d2y2]]
                                        for g in combinations_to_try:
                                            if d1x1 < g[0] < d1x2 and d1y1 < g[1] < d1y2:
                                                if ("VH" in str(strings[i][j]) and "VL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("VL" in str(strings[i][j]) and "VH" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CL" in str(strings[i][j]) and "CH1" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH1" in str(strings[i][j]) and "CL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH2" in str(strings[i][j]) and "CH2" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH3" in str(strings[i][j]) and "CH3" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH4" in str(strings[i][j]) and "CH4" in str(domains_dict.get(domains_keyslist[f])[1])) or ("-H-" == str(strings[i][j]) and "-H-" == str(domains_dict.get(domains_keyslist[f])[1])) or ("X" in str(strings[i][j]) and "X" in str(domains_dict.get(domains_keyslist[f])[1])) or ( str(strings[i][j]) == "C" and  str(domains_dict.get(domains_keyslist[f])[1]) == "C"):
                                                    disulphide_count = 0
                                                    paired_X_domains.append(int(paired_number))
                                                    for y in range(len(disulphides_keyslist)):
                                                            #print("Looking for those disulphides")
                                                        disulphx1 = disulphides_dict.get(disulphides_keyslist[y])[0][0]
                                                        disulphy1 = disulphides_dict.get(disulphides_keyslist[y])[0][1]
                                                        disulphx2 = disulphides_dict.get(disulphides_keyslist[y])[0][2]
                                                        disulphy2 = disulphides_dict.get(disulphides_keyslist[y])[0][3]
                                                            #print("disulphx1",d1x1, disulphx1, d1x2 ,"disulphy1", d1y1, disulphy1 , d1y2,"disulphx2" , d2x1 , disulphx2 , d2x2 ,"disulphy2", d2y1, disulphy2 , d2y2)
                                                            #print("disulphx2",d1x1, disulphx2, d1x2 ,"disulphy2", d1y1, disulphy2 , d1y2,"disulphx1" , d2x1 , disulphx1 , d2x2 ,"disulphy1", d2y1, disulphy1 , d2y2)

                                                        if ((d1x1 <= disulphx2 <= d1x2 and d1y1 <= disulphy2 <= d1y2) and (d2x1 <= disulphx1 <= d2x2 and d2y1 <= disulphy1 <= d2y2)) or ((d1x1 <= disulphx1 <= d1x2 and d1y1 <= disulphy1 <= d1y2) and (d2x1 <= disulphx2 <= d2x2 and d2y1 <= disulphy2 <= d2y2)):
                                                            disulphide_count += 1



                    if len(paired_X_domains) > 1:
                        set_paired = set(paired_X_domains)
                        list_set_paired = list(set_paired)
                    else:
                        list_set_paired = paired_X_domains
                    list_paired_X_domains = list_set_paired
                    if len(list_paired_X_domains) > 0:
                        strings[i][j] = str(domain_name+"("+str(number)+":"+str(list_paired_X_domains)+")")
                        strings[i][j] = re.sub("\[|\]","",strings[i][j])
                        for f in range(len(list_paired_X_domains)):
                            for a in range(len(full_chains)):
                                for b in range(len(full_chains[a])):
                                    if ":" not in strings[a][b] and ("X" in strings[a][b] or strings[a][b] == "C") and "-" not in strings[a][b]  :
                                        #if int(full_chains[a][b]) == int(list_paired_X_domains[f]) and ":" not in strings[a][b]:
                                        paired_domain = strings[a][b]
                                        paired_number =  re.findall("\((.*?)\)", str(paired_domain))
                                        paired_number =  str(re.sub("\[|\'|\]","", str(paired_number)))
                                        paired_name   = re.sub("\((.*?)\)", "", str(paired_domain))
                                        if int(paired_number) == int(list_paired_X_domains[f]):
                                            list_paired_X_domains.append(number)
                                            list_to_add = []
                                            for k in list_paired_X_domains:
                                                if int(k) != int(paired_number) and int(k) not in list_to_add:
                                                    list_to_add.append(int(k))
                                            list_to_add_sorted = list_to_add.sort()
                                            strings[a][b] = str(paired_name+"("+str(paired_number)+":"+str(list_to_add)+")")
                                            strings[a][b] = re.sub("\[|\]","",strings[a][b])
                                            list_paired_X_domains.remove(int(number))

                                                        #if disulphide_count == 0:
                                                        #    strings[i][j] = str(domain_name+"("+str(number)+":"+str(paired_number)+")")
                                                        #    strings[a][b] = str(paired_name+"("+str(paired_number)+":"+str(number)+")")
                                                        #elif disulphide_count > 0:
                                                        #    strings[i][j] = str(domain_name+"("+str(number)+":"+str(paired_number)+")"+"{"+str(disulphide_count)+"}")
                                                        #    strings[a][b] = str(paired_name+"("+str(paired_number)+":"+str(number)+")"+"{"+str(disulphide_count)+"}")
                                                        #paired.append(int(number))
                                                        #paired.append(int(paired_number))

            elif ":" not in str(strings[i][j]) and "-H(" in str(strings[i][j]):
                number =  re.findall("\((.*?)\)", str(strings[i][j]))
                number =  int(re.sub("\[|\'|\]","", str(number)))

                if number not in paired:
                    domain_name = re.sub("\((.*?)\)","",str(strings[i][j]))
                    index = full_chains[i][j]
                    d1x1 = (bonds_dict.get(index)[0][0])
                    d1x2 = (bonds_dict.get(index)[0][2])
                    if d1x1 > d1x2:
                        d1x1 = (bonds_dict.get(index)[0][2])-100
                        d1x2 = (bonds_dict.get(index)[0][0])+100
                    else:
                        d1x1 = (bonds_dict.get(index)[0][0])-100
                        d1x2 = (bonds_dict.get(index)[0][2])+100
                    d1y1 = (bonds_dict.get(index)[0][1])-5
                    d1y2 = (bonds_dict.get(index)[0][3])+5

                    for f in range(len(bonds_keyslist)):
                        if bonds_keyslist[f] != index and "-H-" in bonds_dict.get(bonds_keyslist[f])[1]:
                            for a in range(len(full_chains)):
                                for b in range(len(full_chains[a])):
                                    if int(full_chains[a][b]) == int(bonds_keyslist[f]) and ":" not in strings[a][b]:
                                        paired_domain = strings[a][b]
                                        paired_number =  re.findall("\((.*?)\)", str(paired_domain))
                                        paired_number =  str(re.sub("\[|\'|\]","", str(paired_number)))
                                        paired_name   = re.sub("\((.*?)\)", "", str(paired_domain))
                                        if int(paired_number) not in paired:
                                            d2x1 = (bonds_dict.get(bonds_keyslist[f])[0][0])
                                            d2x2 = (bonds_dict.get(bonds_keyslist[f])[0][2])
                                            if d2x1 > d1x2:
                                                d2x1 = (bonds_dict.get(bonds_keyslist[f])[0][2])-150
                                                d2x2 = (bonds_dict.get(bonds_keyslist[f])[0][0])+150
                                            else:
                                                d2x1 = (bonds_dict.get(bonds_keyslist[f])[0][0])-150
                                                d2x2 = (bonds_dict.get(bonds_keyslist[f])[0][2])+150
                                            d2y1 = (bonds_dict.get(bonds_keyslist[f])[0][1])-5
                                            d2y2 = (bonds_dict.get(bonds_keyslist[f])[0][3])+5
                                            combinations_to_try = [[d2x1,d2y1],[d2x2,d2y1],[d2x2,d2y2],[d2x1,d2y2],[d2x1,(d2y1+d2y2/2)],[d2x2,(d2y1+d2y2/2)]]
                                            for g in combinations_to_try:
                                                if d1x1 <= g[0] <= d1x2 and d1y1 <= g[1] <= d1y2:
                                                    disulphide_count = 0
                                                    for y in range(len(disulphides_keyslist)):
                                                        #print("Looking for those disulphides")
                                                        disulphx1 = disulphides_dict.get(disulphides_keyslist[y])[0][0]
                                                        disulphy1 = disulphides_dict.get(disulphides_keyslist[y])[0][1]
                                                        disulphx2 = disulphides_dict.get(disulphides_keyslist[y])[0][2]
                                                        disulphy2 = disulphides_dict.get(disulphides_keyslist[y])[0][3]
                                                        #print("disulphx1",d1x1, disulphx1, d1x2 ,"disulphy1", d1y1, disulphy1 , d1y2,"disulphx2" , d2x1 , disulphx2 , d2x2 ,"disulphy2", d2y1, disulphy2 , d2y2)
                                                        #print("disulphx2",d1x1, disulphx2, d1x2 ,"disulphy2", d1y1, disulphy2 , d1y2,"disulphx1" , d2x1 , disulphx1 , d2x2 ,"disulphy1", d2y1, disulphy1 , d2y2)

                                                        if ((d1x1 <= disulphx2 <= d1x2 and d1y1<= disulphy2 <= d1y2) and (d2x1 <= disulphx1 <= d2x2 and d2y1<= disulphy1 <= d2y2)) or ((d1x1 <= disulphx1 <= d1x2 and d1y1<= disulphy1 <= d1y2) and (d2x1 <= disulphx2 <= d2x2 and d2y1 <= disulphy2 <= d2y2)):
                                                            disulphide_count += 1

                                                    if disulphide_count == 0:
                                                        strings[i][j] = str("-H"+"("+str(number)+":"+str(paired_number)+")-")
                                                        strings[a][b] = str("-H"+"("+str(paired_number)+":"+str(number)+")-")
                                                    elif disulphide_count > 0:
                                                        strings[i][j] = str("-H"+"("+str(number)+":"+str(paired_number)+")"+"{"+str(disulphide_count)+"}-")
                                                        strings[a][b] = str("-H"+"("+str(paired_number)+":"+str(number)+")"+"{"+str(disulphide_count)+"}-")
                                                    paired.append(int(number))
                                                    paired.append(int(paired_number))


##Find comments on domains and not on domains
    for i in range(len(full_chains)):
        for j in range(len(full_chains[i])):
            #print(domains_dict.get(full_chains[i][j])[0])
            if "-" not in strings[i][j]:
                domain_type = strings[i][j].split("(")[0]
                domain_type = re.sub("\.|\*|\+|\-|\@|\>","", str(domain_type))
                coordinates = domains_dict.get(full_chains[i][j])[0]
                min_max = get_min_max_coordinates(coordinates)
                d1x1 = min_max[0]
                d1x2 = min_max[1]
                d1y1 = min_max[2]
                d1y2 = min_max[3]
                for k in range(len(type_keyslist)):
                    comment = type_dict.get(type_keyslist[k])[1]
                    labelx = type_dict.get(type_keyslist[k])[0][0]
                    labely = type_dict.get(type_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        note = type_dict.get(type_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("[TYPE:"+note+"]")
                        strings[i][j] += noting
                for k in range(len(note_keyslist)):
                    comment =note_dict.get(note_keyslist[k])[1]
                    labelx = note_dict.get(note_keyslist[k])[0][0]
                    labely = note_dict.get(note_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        note = note_dict.get(note_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("[NOTE:"+note+"]")
                        strings[i][j] += noting
                for k in range(len(mod_keyslist)):
                    comment =mod_dict.get(mod_keyslist[k])[1]
                    labelx = mod_dict.get(mod_keyslist[k])[0][0]
                    labely = mod_dict.get(mod_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        note = mod_dict.get(mod_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("[MOD:"+note+"]")
                        strings[i][j] += noting
                for k in range(len(anti_keyslist)):
                    comment =anti_dict.get(anti_keyslist[k])[1]
                    labelx = anti_dict.get(anti_keyslist[k])[0][0]
                    labely = anti_dict.get(anti_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        note = anti_dict.get(anti_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("[ANTI:"+note+"]")
                        strings[i][j] += noting
                for k in range(len(length_keyslist)):
                    comment =length_dict.get(length_keyslist[k])[1]
                    labelx = length_dict.get(length_keyslist[k])[0][0]
                    labely = length_dict.get(length_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        note = length_dict.get(length_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("[LENGTH:"+note+"]")
                        strings[i][j] += noting
            elif "-H" in strings[i][j]:
                domain_type = strings[i][j].split("(")[0]
                domain_type = re.sub("\.|\*|\+|\-|\@|\>","", str(domain_type))
                coordinates = bonds_dict.get(full_chains[i][j])[0]
                min_max = get_min_max_coordinates(coordinates)
                d1x1 = min_max[0]-50
                d1x2 = min_max[1]+50
                d1y1 = min_max[2]
                d1y2 = min_max[3]
                for k in range(len(type_keyslist)):
                    comment = type_dict.get(type_keyslist[k])[1]
                    labelx = type_dict.get(type_keyslist[k])[0][0]
                    labely = type_dict.get(type_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = type_dict.get(type_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("-"+domain+"[TYPE:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(note_keyslist)):
                    comment =note_dict.get(note_keyslist[k])[1]
                    labelx = note_dict.get(note_keyslist[k])[0][0]
                    labely = note_dict.get(note_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = note_dict.get(note_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        domain = re.sub("\(","*(", domain)
                        noting = str("-"+domain+"[NOTE:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(mod_keyslist)):
                    comment =mod_dict.get(mod_keyslist[k])[1]
                    labelx = mod_dict.get(mod_keyslist[k])[0][0]
                    labely = mod_dict.get(mod_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = mod_dict.get(mod_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        domain = re.sub("\(","*(", domain)
                        noting = str("-"+domain+"[MOD:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(anti_keyslist)):
                    comment =anti_dict.get(anti_keyslist[k])[1]
                    labelx = anti_dict.get(anti_keyslist[k])[0][0]
                    labely = anti_dict.get(anti_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = anti_dict.get(anti_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("-"+domain+"[ANTI:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(length_keyslist)):
                    comment =length_dict.get(length_keyslist[k])[1]
                    labelx = length_dict.get(length_keyslist[k])[0][0]
                    labely = length_dict.get(length_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = length_dict.get(length_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("-"+domain+"[LENGTH:"+note+"]-")
                        strings[i][j] = noting
            elif "-L" in strings[i][j]:
                domain_type = "Linker"
                coordinates = bonds_dict.get(full_chains[i][j])[0]
                min_max = get_min_max_coordinates(coordinates)
                d1x1 = min_max[0]-50
                d1x2 = min_max[1]+50
                d1y1 = min_max[2]
                d1y2 = min_max[3]
                for k in range(len(type_keyslist)):
                    comment = type_dict.get(type_keyslist[k])[1]
                    labelx = type_dict.get(type_keyslist[k])[0][0]
                    labely = type_dict.get(type_keyslist[k])[0][1]
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = type_dict.get(type_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("-L[TYPE:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(note_keyslist)):
                    comment =note_dict.get(note_keyslist[k])[1]
                    labelx = note_dict.get(note_keyslist[k])[0][0]
                    labely = note_dict.get(note_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = note_dict.get(note_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        domain = re.sub("\(","*(", domain)
                        noting = str("-L[NOTE:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(mod_keyslist)):
                    comment =mod_dict.get(mod_keyslist[k])[1]
                    labelx = mod_dict.get(mod_keyslist[k])[0][0]
                    labely = mod_dict.get(mod_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = mod_dict.get(mod_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        domain = re.sub("\(","*(", domain)
                        noting = str("-L[MOD:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(anti_keyslist)):
                    comment =anti_dict.get(anti_keyslist[k])[1]
                    labelx = anti_dict.get(anti_keyslist[k])[0][0]
                    labely = anti_dict.get(anti_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        domain = strings[i][j].split("-")[1]
                        note = anti_dict.get(anti_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("-L[ANTI:"+note+"]-")
                        strings[i][j] = noting
                for k in range(len(length_keyslist)):
                    comment =length_dict.get(length_keyslist[k])[1]
                    labelx = length_dict.get(length_keyslist[k])[0][0]
                    labely = length_dict.get(length_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or domain_type in str(comment):
                        #domain = strings[i][j].split("-")[1]
                        note = length_dict.get(length_keyslist[k])[1]
                        note = re.sub(domain_type, "",note)
                        note = re.sub(" NOTE:| TYPE:| MOD:| ANTI:| LENGTH:","", note)
                        noting = str("-L[LENGTH:"+note+"]-")
                        strings[i][j] = noting

##conver lists to expression
    final_string = ""
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            if j == 0 and i == 0:
                final_string += strings[i][j]
            elif j == 0 and i > 0:
                final_string+= str("|"+strings[i][j])
            elif j > 0 :
                final_string += strings[i][j]
    final_string = re.sub("\-\-","-", str(final_string))
    final_string = re.sub("\-L\-\|","|",str(final_string))
    final_string = re.sub("\-\|","|",str(final_string))
    final_string = re.sub("\-$","",str(final_string))


##display string to textbox
    textBox.delete("1.0","end")
    textBox.insert("1.0",str(final_string))


def button_hover(e):
    button["bg"] = "white"
    status_label.config(text="Get Input")
def button_hover_leave(e):
    button["bg"] = "grey"
    status_label.config(text="")
def button_hover_polygon(label):
    status_label.config(text=label)
def button_hover_polygon_leave(e):
    status_label.config(text="")
def delete_all_button(canvas):
    lower_canvas.delete("all")
    global canvas_polygons
    global canvas_labels
    global temp_label
    global TYPE_labels
    global NOTE_labels
    global domain_type
    global domain_mod
    global domain_charge
    global Delete_lock
    global Bond_lock
    global Domain_Primer
    global Domain_Primer_Lock
    global extra_mods
    global all_buttons
    for i in range(len(all_buttons)):
        all_buttons[i].config(fg = "black")
    canvas_polygons = {}
    canvas_labels = {}
    temp_label    = {}
    TYPE_labels   = {}
    NOTE_labels   = {}
    MOD_labels    = {}
    ANTI_labels   = {}
    LENGTH_labels = {}
    Delete_lock   = ""
    domain_type   = ""
    domain_charge = ""
    domain_mod    = ""
    extra_domains = ""
    Bond_lock     = ""
    Domain_Primer = []
    Domain_Primer_Lock=""
    textBox.delete("1.0","end")


def delete_button(canvas):
    '''
    Mouse click deletes object from canvas
    '''
    global Delete_lock
    global Bond_lock
    global InsertDelClickButton
    global Domain_Primer
    global Domain_Primer_Lock
    global extra_mods
    global domain_type
    global domain_charge
    global domain_mod
    global all_buttons
    Bond_lock = ""
    domain_type = ""
    domain_mod = ""
    domain_charge = ""
    extra_mods = ""
    Domain_Primer = []
    Domain_Primer_Lock = ""

    for i in all_buttons:
        i.config(fg = "black")
    if Delete_lock == False:
        InsertDelClickButton.config(fg="red")
        Delete_lock = True
        lower_canvas.config(cursor = "arrow")
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>",mm.delete)
        status_label.config(text="Delete lock on")

    elif Delete_lock == True:
        Delete_lock = False
        InsertDelClickButton.config(fg="black")
        lower_canvas.unbind("<Button-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")

def domain_type_button(letter):
    global Bond_lock
    global Delete_lock
    global Label_lock
    global domain_charge
    global domain_mod
    global extra_mods
    global domain_type
    global Domain_Primer_Lock
    global Domain_Primer
    global bond_buttons
    global comments_buttons
    global domain_buttons
    global delete_buttons
    global specificity_buttons
    for i in bond_buttons:
        i.config(fg = "black")
    for i in comments_buttons:
        i.config(fg = "black")
    for i in delete_buttons:
        i.config(fg = "black")
    #for i in specificity_buttons:
    #    i.config(fg = "black")
    Delete_lock = False
    Bond_lock = ""
    if Domain_Primer_Lock != "":
        lower_canvas.bind("<Button-1>", mm.place_domain)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        lower_canvas.config(cursor = "plus")
    if Domain_Primer_Lock == "":
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
    status_label.config(text="")
    domain_letter = str(letter)

    if letter not in str(domain_type):
        domain_type += letter
        if letter == "a":
            a_button.config(fg="red")
        elif letter == "b":
            b_button.config(fg="red")
        elif letter == "c":
            c_button.config(fg="red")
        elif letter == "d":
            d_button.config(fg="red")
        elif letter == "e":
            e_button.config(fg="red")
        elif letter == "f":
            f_button.config(fg="red")
        elif letter == "g":
            g_button.config(fg="red")
        elif letter == "h":
            h_button.config(fg="red")
        if Domain_Primer != []:
            Domain_Primer[8] = re.sub("a|b|c|d|e|f|g|h","",Domain_Primer[8])
            Domain_Primer[8] = re.sub("\.","."+domain_type,Domain_Primer[8])
        elif Domain_Primer_Lock == "":
            lower_canvas.config(cursor = "arrow")
            lower_canvas.bind("<Button-1>", mm.change_specificity)
            lower_canvas.unbind("<B1-Motion>")
            lower_canvas.unbind("<ButtonRelease-1>")

    elif letter in str(domain_type): #and Domain_Primer_Lock == "":
        domain_type = re.sub(str(letter),"",domain_type)
        print(domain_type)

        Domain_Primer_Lock = ""
        Domain_Primer = []
        for i in domain_buttons:
            i.config(fg="black")


        if letter == "a":
            a_button.config(fg="black")
        elif letter == "b":
            b_button.config(fg="black")
        elif letter == "c":
            c_button.config(fg="black")
        elif letter == "d":
            d_button.config(fg="black")
        elif letter == "e":
            e_button.config(fg="black")
        elif letter == "f":
            f_button.config(fg="black")
        elif letter == "g":
            g_button.config(fg="black")
        elif letter == "h":
            h_button.config(fg="black")
        if domain_type == "" and Domain_Primer_Lock != "":
            lower_canvas.bind("<Button-1>", mm.select)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.release)
            lower_canvas.config(cursor = "arrow")
        elif domain_type != "" and Domain_Primer_Lock == "":
            if (domain_charge != "" or domain_mod !="" or extra_mods != ""):
                lower_canvas.config(cursor = "arrow")
                lower_canvas.bind("<Button-1>", mm.change_modification)
                lower_canvas.unbind("<B1-Motion>")
                lower_canvas.unbind("<ButtonRelease-1>")
            elif (domain_charge == "" and domain_mod =="" and extra_mods == ""):
                lower_canvas.config(cursor = "arrow")
                lower_canvas.bind("<Button-1>", mm.change_specificity)
                lower_canvas.unbind("<B1-Motion>")
                lower_canvas.unbind("<ButtonRelease-1>")
        else:
            lower_canvas.bind("<Button-1>", mm.select)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.release)
            lower_canvas.config(cursor = "arrow")

def extra_mod_button(letter):
    global extra_mods
    global domain_type
    global domain_charge
    global domain_mod
    global Bond_lock
    global Delete_lock
    global Label_lock
    global Domain_Primer
    global Domain_Primer_Lock
    Delete_lock = False
    Bond_lock = ""
    global bond_buttons
    global delete_buttons
    global comments_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")
    if Domain_Primer_Lock != "":
        lower_canvas.bind("<Button-1>", mm.place_domain)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        lower_canvas.config(cursor = "plus")
    elif Domain_Primer_Lock == "":
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
    selected_mod = str(letter)
    if selected_mod not in str(extra_mods):
        extra_mods += selected_mod
        if selected_mod == "!":
            Gly_button.config(fg = "red")
        elif selected_mod == "*":
            Mod_button.config(fg="red")
        if Domain_Primer != [] and (selected_mod == "!" or selected_mod == "*"):
            Domain_Primer[8] = re.sub("\!|\*","",Domain_Primer[8])
            if "V" in Domain_Primer[8]:
                Domain_Primer[8] = re.sub("\.",extra_mods+".",Domain_Primer[8])
            else:
                Domain_Primer[8] += extra_mods

    elif selected_mod in str(extra_mods):
        if selected_mod == "!":
            extra_mods = re.sub("\!","",extra_mods)
            Gly_button.config(fg = "black")
            if Domain_Primer != []:
                Domain_Primer[8] = re.sub("\!","",Domain_Primer[8])
        elif selected_mod == "*":
            extra_mods = re.sub("\*","",extra_mods)
            Mod_button.config(fg="black")
            if Domain_Primer !=[]:
                Domain_Primer[8] = re.sub("\*","",Domain_Primer[8])

    if (domain_charge != "" or domain_mod !="" or extra_mods != "") and domain_type == "" and Domain_Primer_Lock == "":
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<Button-1>", mm.change_modification)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")


def domain_mod_button(letter):
    domain_letter = str(letter)
    global domain_mod
    global domain_charge
    global extra_mods
    global domain_type
    global Bond_lock
    global Delete_lock
    global Label_lock
    global Domain_Primer
    global Domain_Primer_Lock
    Delete_lock = False
    Bond_lock = ""
    global bond_buttons
    global delete_buttons
    global comments_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")

    if Domain_Primer_Lock != "":
        lower_canvas.bind("<Button-1>", mm.place_domain)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        lower_canvas.config(cursor = "plus")
    elif Domain_Primer_Lock == "":
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
    buttons = [KIH_hole, KIH_knob]
    for i in buttons:
        i.config(fg="black")
    if domain_mod != str(letter):
        domain_mod = str(letter)
        if letter == ">" or letter == "@":
            global domain_direction
            domain_direction = "innie"
            if letter == ">":
                KIH_knob.configure(fg="red")
            elif letter == "@":
                KIH_hole.configure(fg="red")

        if Domain_Primer != [] and (letter == "@" or letter == ">"):
            Domain_Primer[5] = re.sub(Domain_Primer[5],letter,Domain_Primer[5])
            re.sub("\@|\>","", Domain_Primer[8])
            if "." in Domain_Primer[8]:
                Domain_Primer[8] = re.sub("\.", letter+".",Domain_Primer[8])
            else:
                Domain_Primer[8] += letter

    elif domain_mod == str(letter) and Domain_Primer_Lock != "":
        domain_mod = ""
        if Domain_Primer != []:
            Domain_Primer[5] = re.sub(".*",letter,Domain_Primer[5])
            re.sub("\@|\>","", Domain_Primer[8])
            if "." in Domain_Primer[8]:
                Domain_Primer[8] = re.sub("\.", letter+".",Domain_Primer[8])
            else:
                Domain_Primer[8] += letter

    elif domain_mod == str(letter):
        domain_mod = ""
        if Domain_Primer != []:
            Domain_Primer[8] = re.sub("\@|\>","", Domain_Primer[8])
    if (domain_charge != "" or domain_mod != "" or extra_mods != "") and domain_type == "" and Domain_Primer_Lock == "":
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<Button-1>", mm.change_modification)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")



def domain_charge_button(charge):
    letter = str(charge)
    global domain_mod
    global Bond_lock
    global Delete_lock
    global Label_lock
    global domain_charge
    global Domain_Primer_Lock
    global Domain_Primer
    Delete_lock = False
    Bond_lock = ""
    if Domain_Primer_Lock != "":
        lower_canvas.bind("<Button-1>", mm.place_domain)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        lower_canvas.config(cursor = "plus")
    elif Domain_Primer_Lock == "":
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
    global bond_buttons
    global delete_buttons
    global comments_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")
    buttons = [Positive_charge, Negative_charge]
    for i in buttons:
        i.config(fg="black")
    if domain_charge != str(letter):
        domain_charge = str(letter)
        if letter == "+":
            Positive_charge.configure(fg="red")
        elif letter == "_":
            Negative_charge.configure(fg="red")
        if Domain_Primer_Lock != "" and (letter == "+" or letter == "_"):
            Domain_Primer[8] = re.sub("\+|\-|\_","",Domain_Primer[8])
            if "V" in Domain_Primer[8]:
                Domain_Primer[8] = re.sub("\.",letter+".",Domain_Primer[8])
            else:
                Domain_Primer[8] += letter

    elif domain_charge == str(letter):
        domain_charge = ""
        if Domain_Primer != []:
            Domain_Primer[8] = re.sub("\+|\_|\-","",Domain_Primer[8])

    if (domain_charge != "" or domain_mod !="" or extra_mods != "") and domain_type == "" and Domain_Primer_Lock == "":
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<Button-1>", mm.change_modification)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")

def bond_drag_button(canvas,name,buttonpress):
    global Bond_lock
    global Delete_lock
    Delete_lock = False
    global delete_buttons
    global domain_buttons
    global comments_buttons
    global mod_buttons
    global specificity_buttons
    global bond_buttons
    for i in delete_buttons:
        i.config(fg="black")
    for i in domain_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")
    for i in specificity_buttons:
        i.config(fg="black")
    for i in mod_buttons:
        i.config(fg="black")
    for i in bond_buttons:
        i.config(fg="black")
    if Bond_lock != buttonpress:
        print("Connector")
        Bond_lock = buttonpress
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.start_bond)
        lower_canvas.bind("<B1-Motion>", mm.drag_bond)
        InsertBondButton.config(fg="red")
        lower_canvas.bind("<ButtonRelease-1>", mm.release_bond)
        lower_canvas.config(cursor = "cross")

    elif str(Bond_lock) == str(buttonpress):
        Bond_lock = ""
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")
        lower_canvas.config(cursor = "arrow")

def Linker_bond_button(canvas,name,buttonpress):
    global Bond_lock
    global Delete_lock
    global delete_buttons
    global domain_buttons
    global comments_buttons
    global mod_buttons
    global specificity_buttons
    global bond_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in domain_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")
    for i in specificity_buttons:
        i.config(fg="black")
    for i in mod_buttons:
        i.config(fg="black")
    Delete_lock = False
    if Bond_lock != buttonpress:
        Bond_lock = buttonpress
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.start_bond)
        lower_canvas.bind("<B1-Motion>", mm.drag_Linker_bond)
        InsertLinkerButton.config(fg="red")
        lower_canvas.bind("<ButtonRelease-1>", mm.release_Linker_bond)
        status_label.config(text="Hinge lock on")
        lower_canvas.config(cursor = "cross")

    elif str(Bond_lock) == str(buttonpress):
        Bond_lock = ""
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")
        lower_canvas.config(cursor = "arrow")


def disulphide_bond_button(canvas,name,buttonpress):
    global Bond_lock
    global Delete_lock
    global delete_buttons
    global domain_buttons
    global comments_buttons
    global mod_buttons
    global specificity_buttons
    global bond_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in domain_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")
    for i in specificity_buttons:
        i.config(fg="black")
    for i in mod_buttons:
        i.config(fg="black")
    Delete_lock = False
    if Bond_lock != buttonpress:
        Bond_lock = buttonpress
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.start_bond)
        lower_canvas.bind("<B1-Motion>", mm.drag_disulphide_bond)
        InsertDBondButton.config(fg="red")
        if name == "-":
            lower_canvas.bind("<ButtonRelease-1>", mm.release_Disulphide_bond)
            status_label.config(text="Disulphide lock on")
        lower_canvas.config(cursor = "cross")

    elif str(Bond_lock) == str(buttonpress):
        Bond_lock = ""
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")
        lower_canvas.config(cursor = "arrow")

def Hinge_bond_button(canvas,name,buttonpress):
    global Bond_lock
    global Delete_lock
    global delete_buttons
    global domain_buttons
    global comments_buttons
    global mod_buttons
    global specificity_buttons
    global bond_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in domain_buttons:
        i.config(fg="black")
    for i in comments_buttons:
        i.config(fg="black")
    for i in specificity_buttons:
        i.config(fg="black")
    for i in mod_buttons:
        i.config(fg="black")
    Delete_lock = False
    if Bond_lock != buttonpress:
        Bond_lock = buttonpress
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.start_bond)
        lower_canvas.bind("<B1-Motion>", mm.drag_Hinge_bond)
        InsertLHingeButton.config(fg="red")
        lower_canvas.bind("<ButtonRelease-1>", mm.release_Hinge_bond)
        status_label.config(text="Hinge lock on")
        lower_canvas.config(cursor = "cross")

    elif str(Bond_lock) == str(buttonpress):
        Bond_lock = ""
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")
        lower_canvas.config(cursor = "arrow")

def labels_button(canvas):
    global Label_lock
    global canvas_labels
    global canvas_polygons
    global temp_label
    if Label_lock == False:
        Label_lock = True
    elif Label_lock == True:
        Label_lock = False

    label_keyslist   =list(canvas_labels.keys())
    polygons_keyslist=list(canvas_polygons.keys())
    #DELETE LABELS AND RE RENDER
    for i in range(len(label_keyslist)):
        del canvas_labels[label_keyslist[i]]
        lower_canvas.delete(label_keyslist[i])
    if Label_lock == True:
        for i in range(len(polygons_keyslist)):
            domain_label = str(canvas_polygons.get(polygons_keyslist[i])[1])
            domain_coordinates= canvas_polygons.get(polygons_keyslist[i])[0]
            if "-" not in str(domain_label):
                domain_name = re.sub("\.|@|>","",domain_label)
                min_max = get_min_max_coordinates(canvas_polygons.get(polygons_keyslist[i])[0])
                labelx = min_max[4]
                labely = min_max[5]
                label  = lower_canvas.create_text([labelx,labely], text = domain_name, tags = "label")
                canvas_labels[label] = [[labelx, labely], domain_name]
                temp_label = {}
            elif domain_label == "-H-":
                domain_name = re.sub("-","", domain_label)
                firstx  = domain_coordinates[0]
                secondx = domain_coordinates[2]
                firsty  = domain_coordinates[1]
                secondy = domain_coordinates[3]
                if firstx < secondx:
                    labelx = (firstx) - 30
                elif firstx > secondx:
                    labelx = (firstx) + 30
                labely = (firsty+secondy)/2
                label  = lower_canvas.create_text([labelx,labely], text = domain_name, tags = "label")
                canvas_labels[label] = [[labelx, labely], domain_name]
                temp_label = {}
def SelectCommentTypeButton(letter):
    global CustomLabelLock
    global domain_mod
    global Bond_lock
    global Delete_lock
    global Label_lock
    global domain_charge
    global domain_mod
    global Domain_Primer
    global Domain_Primer_Lock
    Delete_lock = False
    Bond_lock = ""
    Domain_Primer = []
    Domain_Primer_Lock = ""
    domain_type = ""
    extra_mods = ""
    domain_mod = ""
    domain_charge  = ""
    global all_buttons
    for i in all_buttons:
        i.config(fg="black")
    if CustomLabelLock != letter:
        CustomLabelLock = letter
        if letter == "TYPE":
            TypeLabelButton.config(fg="red")
        elif letter == "NOTE":
            NoteLabelButton.config(fg="red")
        elif letter == "MOD":
            ModLabelButton.config(fg="red")
        elif letter == "ANTI":
            AntiLabelButton.config(fg="red")
        elif letter == "LENGTH":
            LengthLabelButton.config(fg="red")
        CommentLabelButton_function(lower_canvas)
    elif CustomLabelLock == letter:
        CustomLabelLock = ""
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")
        lower_canvas.config(cursor = "arrow")

def CommentLabelButton_function(canvas):
    global CustomLabelLock
    global Bond_lock
    global Delete_lock
    global Label_lock
    global Domain_Primer_Lock
    global Domain_Primer
    Delete_lock = False
    Bond_lock = ""
    Domain_Primer_Lock = ""
    Domain_Primer = []
    global delete_buttons
    global domain_buttons
    global mod_buttons
    global specificity_buttons
    global bond_buttons
    for i in bond_buttons:
        i.config(fg="black")
    for i in delete_buttons:
        i.config(fg="black")
    for i in domain_buttons:
        i.config(fg="black")
    for i in specificity_buttons:
        i.config(fg="black")
    for i in mod_buttons:
        i.config(fg="black")

    if CustomLabelLock == "TYPE":
        CustomLabelButton.configure(fg="red")
        entry=CustomLabelEntry.get("1.0","end-1c")
        lower_canvas.bind("<Button-1>", mm.place_type_label)
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        status_label.config(text=entry)
        lower_canvas.config(cursor = "plus")
    elif CustomLabelLock == "NOTE":
        CustomLabelButton.configure(fg="red")
        entry=CustomLabelEntry.get("1.0","end-1c")
        lower_canvas.bind("<Button-1>", mm.place_note_label)
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        status_label.config(text=entry)
        lower_canvas.config(cursor = "plus")
    elif CustomLabelLock == "MOD":
        CustomLabelButton.configure(fg="red")
        entry=CustomLabelEntry.get("1.0","end-1c")
        lower_canvas.bind("<Button-1>", mm.place_mod_label)
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        status_label.config(text=entry)
        lower_canvas.config(cursor = "plus")
    elif CustomLabelLock == "ANTI":
        CustomLabelButton.configure(fg="red")
        entry=CustomLabelEntry.get("1.0","end-1c")
        lower_canvas.bind("<Button-1>", mm.place_anti_label)
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        status_label.config(text=entry)
        lower_canvas.config(cursor = "plus")
    elif CustomLabelLock == "LENGTH":
        CustomLabelButton.configure(fg="red")
        entry=CustomLabelEntry.get("1.0","end-1c")
        lower_canvas.bind("<Button-1>", mm.place_length_label)
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        status_label.config(text=entry)
        lower_canvas.config(cursor = "plus")

    elif CustomLabelLock == "":
        CustomLabelButton.configure(fg="black")
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")
        lower_canvas.config(cursor = "arrow")

def raise_error(canvas,message):
    width = lower_canvas.winfo_width()
    height = lower_canvas.winfo_height()
    lower_canvas.delete("all")
    lower_canvas.create_text((width/2),(height/3), text = message, fill = "red")
    raise IndexError
def items_selected(e):
    '''
    Render item selected in library
    '''
    textBox.delete("1.0","end")
    global Bond_lock
    global Delete_lock
    global antibodyformats
    global formats_keyslist
    global all_buttons
    for i in all_buttons:
        i.config(fg="black")
    Delete_lock = False
    Bond_lock = ""
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "arrow")
    status_label.config(text="")
    i=Library.curselection()
    index = i[0]
    entry = antibodyformats.get(formats_keyslist[index])
    #print(i)
    #print(entry)
    textBox.insert("1.0",str(entry))
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates,lower_canvas,True)

class MouseMover():
    startcoordinates = []
    newcoordinates = []
    label = ""
    def __init__(self):
        self.item = 0; self.previous = (0, 0)
        ###Click and drag item on canvas####
    def select(self, event):
        global canvas_labels
        global temp_label
        global canvas_polygons
        self.startcoordinates = []
        self.newcoordinates =[]
        widget = event.widget
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc, halo = 5)[0]        # ID for closest
        self.previous = (xc, yc)
        #check you haven't selected text and then change self.item
        global canvas_labels
        global temp_label
        global TYPE_labels
        global NOTE_labels
        global ANTI_labels
        global MOD_labels
        global LENGTH_labels
        label_keyslist = list(canvas_labels.keys())
        TYPE_keyslist = list(TYPE_labels.keys())
        NOTE_keyslist = list(NOTE_labels.keys())
        MOD_keyslist = list(MOD_labels.keys())
        ANTI_keyslist = list(ANTI_labels.keys())
        LENGTH_keyslist = list(LENGTH_labels.keys())
        if self.item in label_keyslist:
            polygons_keyslist = list(canvas_polygons.keys())
            for i in range(len(polygons_keyslist)):
                if "-" not in str(canvas_polygons.get(polygons_keyslist[i])[1]):
                    domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
                    min_max = get_min_max_coordinates(domain_coordinates)
                    x1 = min_max[0]
                    x2 = min_max[1]
                    y1 = min_max[2]
                    y2 = min_max[3]
                    if x1< xc <x2 and y1 < yc < y2:
                        self.item = polygons_keyslist[i]

        elif self.item in TYPE_keyslist:
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            return(startcoordinates)
        elif self.item in NOTE_keyslist:
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            return(startcoordinates)
        elif self.item in ANTI_keyslist:
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            return(startcoordinates)
        elif self.item in MOD_keyslist:
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            return(startcoordinates)
        elif self.item in LENGTH_keyslist:
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            return(startcoordinates)

        if "-" not in str(canvas_polygons.get(self.item)[1]):
            domain_coordinates = (canvas_polygons.get(self.item)[0])
            xs= []
            ys= []
            for f in range(len(domain_coordinates)):
                if f%2==0:
                    xs.append(domain_coordinates[f])
                elif f %2 !=0:
                    ys.append(domain_coordinates[f])
            x1 = min(xs)
            x2 = max(xs)
            y1 = min(ys)
            y2 = max(ys)

            if x1< xc <x2 and y1 < yc < y2:
                self.startcoordinates = [xc, yc]
                lower_canvas.config(cursor = "fleur")
                #carry label with domain

                for i in range(len(label_keyslist)):
                    label_text = canvas_labels.get(label_keyslist[i])[1]
                    label_location = canvas_labels.get(label_keyslist[i])[0]
                    labelx = label_location[0]
                    labely = label_location[1]
                    label_text_test = re.sub("\.","",canvas_polygons.get(self.item)[1])
                    if x1<= labelx <=x2 and y1 <= labely <= y2 and label_text_test==label_text:
                        del canvas_labels[label_keyslist[i]]
                        lower_canvas.delete(label_keyslist[i])
                        temp_label = {}
                        temp_label[label_keyslist[i]] = [[labelx,labely], label_text]

                #print("STARTING COORDINATES1", self.startcoordinates)
                return(startcoordinates)
            else:
                self.item = 0
        elif "-" in str(canvas_polygons.get(self.item)[1]):
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            #print("STARTING COORDINATES1", self.startcoordinates)
            return(startcoordinates)

    def drag(self, event):
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        lower_canvas.move(self.item, xc-self.previous[0], yc-self.previous[1])
        self.previous = (xc, yc)
        coordinates = canvas_polygons.get(self.item)
        self.newcoordinates = [xc,yc]
        return(newcoordinates)
    def release(self, event):
        lower_canvas.config(cursor = "arrow")
        global Label_lock
        if self.newcoordinates==[]:
            self.newcoordinates.append(self.startcoordinates[0])
            self.newcoordinates.append(self.startcoordinates[1])
        canvas_keyslist = list(canvas_polygons.keys())
        TYPE_keyslist = list(TYPE_labels.keys())
        if self.item in canvas_keyslist:
            x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
            y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
            diffx = x2-x1
            diffy = y2-y1
            coordinates    = canvas_polygons.get(self.item)[0]
            name           = canvas_polygons.get(self.item)[1]
            new_coordinates= []
            for i in range(len(coordinates)):
                if i%2 ==0:
                    new_coordinates.append((coordinates[i]+diffx))
                elif i%2!=0:
                    new_coordinates.append((coordinates[i]+diffy))
            canvas_polygons[self.item]=[new_coordinates, name]
            print(canvas_polygons.get(self.item))

            if Label_lock == True:
                temp_label_key = list(temp_label.keys())
                if len(temp_label_key) >0:
                    print(temp_label.get(temp_label_key[0])[0])
                    label_name = temp_label.get(temp_label_key[0])[1]
                #label_d1 = temp_label.get(temp_label_key[0])[0]
                firstx  = new_coordinates[2]
                secondx = new_coordinates[4]
                thirdx = new_coordinates[8]
                xs= []
                ys= []
                for f in range(len(new_coordinates)):
                    if f%2==0:
                        xs.append(new_coordinates[f])
                    elif f %2 !=0:
                        ys.append(new_coordinates[f])
                x1 = min(xs)
                x2 = max(xs)
                y1 = min(ys)
                y2 = max(ys)
                #labely2 = (new_coordinates[1]+new_coordinates[5])/2
                labelx2 = (x1+x2)/2
                labely2 = (y1+y2)/2
                label  = lower_canvas.create_text([labelx2,labely2], text = label_name, tags = "label")
                canvas_labels[label] = [[labelx2, labely2], label_name]
                del temp_label[temp_label_key[0]]

        elif self.item in TYPE_keyslist:
            x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
            y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
            diffx = x2-x1
            diffy = y2-y1
            coordinates    = TYPE_labels.get(self.item)[0]
            name           = TYPE_labels.get(self.item)[1]
            new_coordinates= []
            for i in range(len(coordinates)):
                if i%2 ==0:
                    new_coordinates.append((coordinates[i]+diffx))
                elif i%2!=0:
                    new_coordinates.append((coordinates[i]+diffy))
            TYPE_labels[self.item]=[new_coordinates, name]
            print(TYPE_labels.get(self.item))
        print(canvas_polygons)
        print(canvas_labels)
        print(TYPE_labels)
        startcoordinates = []
        newcoordinates = []
    ####Delete selected item on canvas###
    def delete(self,event):
        global deleted_polygons
        global canvas_polygons
        global canvas_labels
        global TYPE_labels
        global NOTE_labels
        global ANTI_labels
        global MOD_labels
        global LENGTH_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc)[0]        # ID for closest
        polygons_keyslist = list(canvas_polygons.keys())
        labels_keyslist = list(canvas_labels.keys())
        TYPE_keyslist = list(TYPE_labels.keys())
        NOTE_keyslist = list(NOTE_labels.keys())
        MOD_keyslist = list(MOD_labels.keys())
        ANTI_keyslist = list(ANTI_labels.keys())
        LENGTH_keyslist = list(LENGTH_labels.keys())
        if self.item in TYPE_keyslist:
            lower_canvas.delete(self.item)
            del TYPE_labels[self.item]
        elif self.item in NOTE_keyslist:
            lower_canvas.delete(self.item)
            del NOTE_labels[self.item]
        elif self.item in ANTI_keyslist:
            lower_canvas.delete(self.item)
            del ANTI_labels[self.item]
        elif self.item in MOD_keyslist:
            lower_canvas.delete(self.item)
            del MOD_labels[self.item]
        elif self.item in LENGTH_keyslist:
            lower_canvas.delete(self.item)
            del LENGTH_labels[self.item]
        elif self.item in labels_keyslist:
            for i in range(len(polygons_keyslist)):
                if "-" not in str(canvas_polygons.get(polygons_keyslist[i])[1]):
                    domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
                    min_max = get_min_max_coordinates(domain_coordinates)
                    x1 = min_max[0]
                    x2 = min_max[1]
                    y1 = min_max[2]
                    y2 = min_max[3]
                    if x1< xc <x2 and y1 < yc < y2:
                        self.item = polygons_keyslist[i]

            domain_coordinates = (canvas_polygons.get(self.item)[0])
            Domain_name = (canvas_polygons.get(self.item)[1])
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            for i in range(len(labels_keyslist)):
                labelx = canvas_labels.get(labels_keyslist[i])[0][0]
                labely = canvas_labels.get(labels_keyslist[i])[0][1]
                label = canvas_labels.get(labels_keyslist[i])[1]
                Domain_name_to_compare = re.sub("\.|\@|\>","",Domain_name)
                Domain_name_to_compare = re.sub("\_","-",Domain_name_to_compare)
                if x1< labelx <x2 and y1 < labely < y2 and str(Domain_name_to_compare) == str(label):
                    lower_canvas.delete(labels_keyslist[i])
                    del canvas_labels[labels_keyslist[i]]
            lower_canvas.delete(self.item)
            del canvas_polygons[self.item]
            deleted_polygons = {self.item:[domain_coordinates,Domain_name]}
        else:
            domain_coordinates = (canvas_polygons.get(self.item)[0])
            Domain_name = (canvas_polygons.get(self.item)[1])
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            for i in range(len(labels_keyslist)):
                labelx = canvas_labels.get(labels_keyslist[i])[0][0]
                labely = canvas_labels.get(labels_keyslist[i])[0][1]
                label = canvas_labels.get(labels_keyslist[i])[1]
                Domain_name_to_compare = re.sub("\.|\@|\>","",Domain_name)
                Domain_name_to_compare = re.sub("\_","-",Domain_name_to_compare)
                if x1< labelx <x2 and y1 < labely < y2 and str(Domain_name_to_compare) == str(label):
                    lower_canvas.delete(labels_keyslist[i])
                    del canvas_labels[labels_keyslist[i]]
            lower_canvas.delete(self.item)
            del canvas_polygons[self.item]
            deleted_polygons = {self.item:[domain_coordinates,Domain_name]}
    ###Click item to reverse orientation###
    def change_specificity(self,event):
        widget = event.widget                       # Get handle to canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc,halo = 5, start="domain")[0]
        global canvas_polygons
        global domain_type
        global canvas_labels
        global temp_label
        global specificity_colours
        label_keyslist = list(canvas_labels.keys())
        polygons_keyslist = list(canvas_polygons.keys())
        if domain_type != "":
            if self.item not in polygons_keyslist:
                if self.item in label_keyslist:
                    polygons_keyslist = list(canvas_polygons.keys())
                    for i in range(len(polygons_keyslist)):
                        if "-" not in str(canvas_polygons.get(polygons_keyslist[i])[1]):
                            domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
                            min_max = get_min_max_coordinates(domain_coordinates)
                            x1 = min_max[0]
                            x2 = min_max[1]
                            y1 = min_max[2]
                            y2 = min_max[3]
                            if x1< xc <x2 and y1 < yc < y2:
                                self.item = polygons_keyslist[i]
            domain_name = canvas_polygons.get(self.item)[1]
            if "V" in domain_name:
                new_domain_name = re.sub("a|b|c|d|e|f|g|h","",domain_name)
                new_domain_name = re.sub("\.","."+domain_type,new_domain_name)
                new_display_name= re.sub("\.|@|>|nano","", new_domain_name)

                if "a" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[0], specificity_colours[1]
                elif "b" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[2], specificity_colours[3]
                elif "c" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[4], specificity_colours[5]
                elif "d" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[6], specificity_colours[7]
                elif "e" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[8], specificity_colours[9]
                elif "f" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[10], specificity_colours[11]
                elif "g" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[12], specificity_colours[13]
                elif "h" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[14], specificity_colours[15]
                elif "X" in str(new_domain_name):
                    heavy_colour, light_colour = specificity_colours[18],specificity_colours[18]
                else:
                    heavy_colour, light_colour = generic_heavy_colour, generic_light_colour
                if "H" in new_domain_name:
                    lower_canvas.itemconfig(self.item, fill = heavy_colour)
                elif "L" in new_domain_name:
                    lower_canvas.itemconfig(self.item, fill = light_colour)
                canvas_polygons[self.item][1] = new_domain_name
                ###change display label
                domain_coordinates = (canvas_polygons.get(self.item)[0])
                min_max = get_min_max_coordinates(domain_coordinates)
                x1 = min_max[0]
                x2 = min_max[1]
                y1 = min_max[2]
                y2 = min_max[3]
                for i in range(len(label_keyslist)):
                    labelx = canvas_labels.get(label_keyslist[i])[0][0]
                    labely = canvas_labels.get(label_keyslist[i])[0][1]
                    if x1 <= labelx <=x2 and y1 <= labely <= y2:
                        del canvas_labels[label_keyslist[i]]
                        lower_canvas.delete(label_keyslist[i])
                        label  = lower_canvas.create_text([labelx,labely], text = new_display_name, tags = "label")
                        canvas_labels[label] = [[labelx,labely], new_display_name]
                        temp_label[label] = [[labelx,labely], new_domain_name]

    def change_modification(self,event):
        widget = event.widget                       # Get handle to canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc,halo = 5, start="domain")[0]
        global canvas_polygons
        global domain_type
        global domain_mod
        global extra_mods
        global domain_charge
        global canvas_labels
        global temp_label
        label_keyslist = list(canvas_labels.keys())
        polygons_keyslist = list(canvas_polygons.keys())
        if domain_mod!= "" or domain_charge != "" or extra_mods != "":
            if self.item not in polygons_keyslist:
                if self.item in label_keyslist:
                    polygons_keyslist = list(canvas_polygons.keys())
                    for i in range(len(polygons_keyslist)):
                        if "-" not in str(canvas_polygons.get(polygons_keyslist[i])[1]):
                            domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
                            min_max = get_min_max_coordinates(domain_coordinates)
                            x1 = min_max[0]
                            x2 = min_max[1]
                            y1 = min_max[2]
                            y2 = min_max[3]
                            if x1< xc <x2 and y1 < yc < y2:
                                self.item = polygons_keyslist[i]
            domain_name = canvas_polygons.get(self.item)[1]
            new_domain_name = re.sub("\@|\>|\+|\-|\_|\!|\*","",domain_name)

            domain_charge = re.sub("\_","-",domain_charge)
            domain_charge_to_add = domain_charge
            if domain_charge in str(domain_name):
                if "+" in domain_charge:
                    domain_charge_to_add = re.sub("\+","",domain_charge_to_add)
                elif "-" in domain_charge:
                    domain_charge_to_add = re.sub("\-","",domain_charge_to_add)

            domain_mod_to_add = domain_mod
            if domain_mod in str(domain_name):
                domain_mod_to_add = re.sub(domain_mod,"",domain_mod_to_add)

            extra_mods_to_add = extra_mods
            if extra_mods in str(domain_name):
                if "*" in extra_mods:
                    extra_mods_to_add = re.sub("\*","",extra_mods_to_add)
                elif "!" in extra_mods:
                    extra_mods_to_add = re.sub("\!","",extra_mods_to_add)

            if "V" in domain_name:
                new_domain_name = re.sub("\.",extra_mods_to_add+domain_mod_to_add+domain_charge_to_add+".",new_domain_name)
                new_display_name= re.sub("\.|nano","", new_domain_name)
            else:
                new_domain_name +=extra_mods_to_add+domain_mod_to_add+domain_charge_to_add
                new_display_name=  re.sub("\.|nano","", new_domain_name)
            canvas_polygons[self.item][1] = new_domain_name
                ###change display label
            domain_coordinates = (canvas_polygons.get(self.item)[0])
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            for i in range(len(label_keyslist)):
                labelx = canvas_labels.get(label_keyslist[i])[0][0]
                labely = canvas_labels.get(label_keyslist[i])[0][1]
                if x1 <= labelx <=x2 and y1 <= labely <= y2:
                    del canvas_labels[label_keyslist[i]]
                    lower_canvas.delete(label_keyslist[i])
                    label  = lower_canvas.create_text([labelx,labely], text = new_display_name, tags = "label")
                    canvas_labels[label] = [[labelx,labely], new_display_name]
                    temp_label[label] = [[labelx,labely], new_domain_name]

    def change_orientation(self,event):
        self.startcoordinates = []
        self.newcoordinates =[]
        widget = event.widget                       # Get handle to canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc,halo = 5, start="domain")[0]
        global canvas_labels
        global temp_label
        global specificity_colours
        label_keyslist = list(canvas_labels.keys())
        if self.item in label_keyslist:
            polygons_keyslist = list(canvas_polygons.keys())
            for i in range(len(polygons_keyslist)):
                if "-" not in str(canvas_polygons.get(polygons_keyslist[i])[1]):
                    domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
                    xs= []
                    ys= []
                    for f in range(len(domain_coordinates)):
                        if f%2==0:
                            xs.append(domain_coordinates[f])
                        elif f %2 !=0:
                            ys.append(domain_coordinates[f])
                    x1 = min(xs)
                    x2 = max(xs)
                    y1 = min(ys)
                    y2 = max(ys)
                    if x1< xc <x2 and y1 < yc < y2:
                        self.item = polygons_keyslist[i]

        coordinates    = canvas_polygons.get(self.item)[0]
        domain_name    = canvas_polygons.get(self.item)[1]
        ##reverse coordinates
        x1 = (canvas_polygons.get(self.item)[0][4])
        x2 = (canvas_polygons.get(self.item)[0][8])
        if x1 > x2:
            x1 = (canvas_polygons.get(self.item)[0][8])
            x2 = (canvas_polygons.get(self.item)[0][4])
        y1 = (canvas_polygons.get(self.item)[0][1])
        y2 = (canvas_polygons.get(self.item)[0][5])
        if x1< xc <x2 and y1 < yc < y2:
            new_coordinates = []
            if "X" not in domain_name:
                for i in range(len(coordinates)):
                    #print(coordinates[i])
                    if i%2 ==0:
                        if coordinates[i] == x1:
                            coordinates[i] = x2
                        elif coordinates[i] == x2:
                            coordinates[i] = x1
                        elif x1 < coordinates[i] and x2 < coordinates[i]:
                            coordinates[i] = x1-20
                        elif x1 > coordinates[i] and x2 > coordinates[i]:
                            coordinates[i] = x2+20
                        new_coordinates.append((coordinates[i]))
                    elif i%2 != 0:
                        new_coordinates.append(coordinates[i])
            elif "X" in domain_name:
                new_coordinates = coordinates
            lower_canvas.delete(self.item)
            del canvas_polygons[self.item]
            new_domain_name = domain_name

            if "a" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[0], specificity_colours[1]
            elif "b" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[2], specificity_colours[3]
            elif "c" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[4], specificity_colours[5]
            elif "d" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[6], specificity_colours[7]
            elif "e" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[8], specificity_colours[9]
            elif "f" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[10], specificity_colours[11]
            elif "g" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[12], specificity_colours[13]
            elif "h" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[14], specificity_colours[15]
            elif "X" in str(new_domain_name):
                heavy_colour, light_colour = specificity_colours[18],specificity_colours[18]
            else:
                heavy_colour, light_colour = generic_heavy_colour, generic_light_colour
            if "VH" in domain_name or "CH" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=heavy_colour, width=2, tags="domain")
            elif "VL" in domain_name or "CL" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=light_colour, width=2, tags="domain")
            elif "X" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=X_colour, width=2, tags="domain")
            elif "C" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=C_colour, width=2, tags="domain")


            if Label_lock == True:
                label_keyslist = list(canvas_labels.keys())
                for i in range(len(label_keyslist)):
                    label_text = canvas_labels.get(label_keyslist[i])[1]
                    label_location = canvas_labels.get(label_keyslist[i])[0]
                    labelx = label_location[0]
                    labely = label_location[1]
                    if x1< labelx <x2 and y1 < labely < y2:
                        del canvas_labels[label_keyslist[i]]
                        lower_canvas.delete(label_keyslist[i])
                        #temp_label[label_keyslist[i]] = [[labelx,labely], label_text]
                        canvas_polygons[domain] = [new_coordinates, domain_name]
                        label = lower_canvas.create_text([labelx,labely], text = label_text, tags = "label")
                        canvas_labels[label] = [[labelx,labely], label_text]
            print(new_coordinates, coordinates)
            print(canvas_polygons)

        else:
            self.item = 0
    ###Place Domains ####
    def place_domain(self,event):
        global Domain_Primer
        global Label_lock
        global specificity_colours
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        startx = xc
        starty = yc-40
        Domain_Primer[5] = re.sub("\+|\-|\_","",Domain_Primer[5])
        if Domain_Primer[5] == "@" or Domain_Primer[5] == ">":
            global domain_direction
            Domain_Primer[3] = "innie"
        All_positions_and_chains = {}
        print(Domain_Primer)
        domaincoordinates = domainmaker(All_positions_and_chains,startx,starty,Domain_Primer[0],Domain_Primer[1],Domain_Primer[2],Domain_Primer[3],Domain_Primer[4],Domain_Primer[5],Domain_Primer[6],Domain_Primer[7])
        domain_name = re.sub("nano","",Domain_Primer[8])
        if "a" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[0], specificity_colours[1]
        elif "b" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[2], specificity_colours[3]
        elif "c" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[4], specificity_colours[5]
        elif "d" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[6], specificity_colours[7]
        elif "e" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[8], specificity_colours[9]
        elif "f" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[10], specificity_colours[11]
        elif "g" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[12], specificity_colours[13]
        elif "h" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[14], specificity_colours[15]
        elif "X" in str(domain_name):
            heavy_colour, light_colour = specificity_colours[18],specificity_colours[18]
        elif (domain_name) == "C":
            heavy_colour, light_colour = specificity_colours[19],specificity_colours[19]
        else:
            heavy_colour, light_colour = generic_heavy_colour, generic_light_colour
        if Domain_Primer[10] == True:
            domain = lower_canvas.create_polygon(domaincoordinates[0], outline='#000000',fill=heavy_colour, width=2, tags="domain")
        elif Domain_Primer[9] == True:
            domain = lower_canvas.create_polygon(domaincoordinates[0], outline='#000000',fill=light_colour, width=2, tags="domain")
        canvas_polygons[domain] = [domaincoordinates[0], domain_name]
        if Label_lock == True:
            domain_name = re.sub("\.|@|>","",domain_name)
            label  = lower_canvas.create_text(domaincoordinates[3], text = str(domain_name), tags = "label")
            canvas_labels[label] = [domaincoordinates[3], domain_name]
        global domain_mod
        #domain_mod = ""
        global domain_direction
        domain_direction = "constant"
    def place_domain_release(self,event):
        lower_canvas.config(cursor = "plus")

    def place_type_label(self,event):
        global TYPE_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "TYPE_labels")
            TYPE_labels[label] = [[xc,yc], entry]
    def place_note_label(self,event):
        global NOTE_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "NOTE_labels")
            NOTE_labels[label] = [[xc,yc], entry]
    def place_mod_label(self,event):
        global MOD_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "MOD_labels")
            MOD_labels[label] = [[xc,yc], entry]
    def place_anti_label(self,event):
        global ANTI_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "ANTI_labels")
            ANTI_labels[label] = [[xc,yc], entry]
    def place_length_label(self,event):
        global LENGTH_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "LENGTH_labels")
            LENGTH_labels[label] = [[xc,yc], entry]
    ###Draw and drag bonds###
    def start_bond(self,event):
        lower_canvas.delete("draggable_line")
        # Convert screen coordinates to canvas coordinates
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.startcoordinates = [xc, yc]
        print(self.startcoordinates)
        return(self.startcoordinates)
    def drag_bond(self,event):
        lower_canvas.delete("draggable_line")
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#000000', width=2,tags = "draggable_line")
        self.newcoordinates = [x2,y2]
        return(newcoordinates)
    def release_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[1], width=2, tags=("bonds","connector"))
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates   = []
    def release_Hinge_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-H-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[2], width=2, tags=("bonds","hinge"))
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordiantes   = []
    def release_Linker_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-L-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[3], width=2, tags=("bonds","linker"))
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates =[]
    def drag_disulphide_bond(self,event):
        lower_canvas.delete("draggable_line")
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[0], width=2,tags = "draggable_line")
        self.newcoordinates = [x2,y2]
    def drag_Hinge_bond(self,event):
        lower_canvas.delete("draggable_line")
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[2], width=2,tags = "draggable_line")
        self.newcoordinates = [x2,y2]
    def drag_Linker_bond(self,event):
        lower_canvas.delete("draggable_line")
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[3], width=2,tags = "draggable_line")
        self.newcoordinates = [x2,y2]

    def release_Disulphide_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-disulphide-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[0], width=2, tags="disulphide")
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates =[]










################Domain Drawer######################
def domainmaker(All_positions_and_chains,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H):
    '''
    Most important function, calculates coordinates needed to draw domain polygons
    '''
    ##First check this domain is not overlapping another that has already been drawn###
    ##Only allow overlap if both domains are "X"
    All_positions_and_chains_list = list(All_positions_and_chains.keys())
    for i in range(len(All_positions_and_chains_list)):
        coordinates = All_positions_and_chains.get(All_positions_and_chains_list[i])[0]
        if len(coordinates) > 2:
            x1 = get_min_max_coordinates(coordinates)[0]
            x2 = get_min_max_coordinates(coordinates)[1]
            y1 = get_min_max_coordinates(coordinates)[2]
            y2 = get_min_max_coordinates(coordinates)[3]
            if startx == coordinates[0] and starty == y1:
                if mod == "X":
                    starty -=80
                elif righthanded == False and direction == "innie":
                    startx -=120
                elif righthanded == True and direction == "innie":
                    startx +=120




    if V == False and  direction == "constant" and mod == "" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx - 20
        elif righthanded == True:
            secondx = firstx + 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx + 20
        elif righthanded == True:
            fourthx     = thirdx - 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx + 20
        elif righthanded == True:
            fifthx     = fourthx - 20
        fifthy         = fourthy
        if slant and righthanded == False:
            sixthx     = fifthx - 40
        elif slant and righthanded == True:
            sixthx     = fifthx + 40
        else:
            sixthx     = fifthx
        sixthy     = fifthy-75
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy]
    elif V == True and direction == "innie" and mod != ">" and mod !="@":
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx - 20
        elif righthanded == True:
            secondx = firstx + 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx + 20
        elif righthanded == True:
            fourthx     = thirdx - 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx + 20
        elif righthanded == True:
            fifthx     = fourthx - 20
        fifthy         = fourthy
        if slant and righthanded == False:
            sixthx     = fifthx - 30
        elif slant and righthanded == True:
            sixthx     = fifthx + 30
        else:
            sixthx     = fifthx
        sixthy     = fifthy-55
        if righthanded == False:
            seventhx = sixthx - 20
        elif righthanded == True:
            seventhx = sixthx + 20
        seventhy    = sixthy
        coordinates = [firstx,firsty,secondx,secondy,thirdx,thirdy,fourthx,fourthy,fifthx,fifthy,sixthx,sixthy,seventhx,seventhy]
    elif V == True and direction == "outie" and mod != ">" and mod !="@":
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx + 20
        elif righthanded == True:
            secondx = firstx - 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx - 20
        elif righthanded == True:
            fourthx     = thirdx + 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx - 20
        elif righthanded == True:
            fifthx     = fourthx + 20
        fifthy         = fourthy
        if slant and righthanded == False:
            sixthx     = fifthx - 30
        elif slant and righthanded == True:
            sixthx     = fifthx + 30
        else:
            sixthx     = fifthx
        sixthy     = fifthy-55
        if righthanded == False:
            seventhx = sixthx + 20
        elif righthanded == True:
            seventhx = sixthx - 20
        seventhy    = sixthy
        coordinates = [firstx,firsty,secondx,secondy,thirdx,thirdy,fourthx,fourthy,fifthx,fifthy,sixthx,sixthy,seventhx,seventhy]
    elif V == True and direction == "Single_Fv_Chain" and  mod != ">" and mod !="@":
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx + 20
        elif righthanded == True:
            secondx = firstx - 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx - 20
        elif righthanded == True:
            fourthx     = thirdx + 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx - 20
        elif righthanded == True:
            fifthx     = fourthx + 20
        fifthy         = fourthy
        if slant and righthanded == False:
            sixthx     = fifthx - 40
        elif slant and righthanded == True:
            sixthx     = fifthx + 40
        else:
            sixthx     = fifthx
        sixthy       = firsty
        if slant == False:
            seventhx    = firstx
        elif slant == True and righthanded == False:
            seventhx    = firstx+10
        elif slant == True and righthanded == True:
            seventhx    = firstx-10
        seventhy    = firsty+20

        eighthx     = secondx
        eighthy     = secondy
        coordinates = [secondx,secondy,secondx,secondy,thirdx,thirdy,fourthx,fourthy,fifthx,fifthy,sixthx,sixthy,seventhx,seventhy]
    elif V == True  and mod == "@" and direction == "outie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx - 20
        elif righthanded == True:
            secondx = firstx + 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx + 20
        elif righthanded == True:
            fourthx     = thirdx - 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx + 20
        elif righthanded == True:
            fifthx     = fourthx - 20
        fifthy         = fourthy
        if slant== True and righthanded == True:
            sixthx     = fourthx+15
        elif slant== True and righthanded == False:
            sixthx     = fourthx-15
        else:
            sixthx     = fourthx
        sixthy     = fifthy-25
        if righthanded == False and slant == True:
            eighthx = firstx+10
        elif righthanded==True and slant == True:
            eighthx = firstx-10
        else:
            eighthx = firstx
        eighthy = firsty+20
        if righthanded==False:
            seventhx=eighthx + 20
        elif righthanded==True:
            seventhx= eighthx - 20
        seventhy= eighthy
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy, eighthx, eighthy]
    elif V == True  and mod == "@" and direction == "innie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx + 20
        elif righthanded == True:
            secondx = firstx - 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx - 20
        elif righthanded == True:
            fourthx     = thirdx + 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx - 20
        elif righthanded == True:
            fifthx     = fourthx + 20
        fifthy         = fourthy
        if slant == True:
            sixthx     = fifthx
        elif slant == False and righthanded == False:
            sixthx     = fifthx + 15
        elif slant == False and righthanded == True:
            sixthx     = fifthx - 15
        sixthy     = fifthy-25
        if righthanded == False and slant == True:
            eighthx = firstx+10
        elif righthanded==True and slant == True:
            eighthx = firstx-10
        else:
            eighthx = firstx
        eighthy = firsty+20
        if righthanded==False:
            seventhx=eighthx - 20
        elif righthanded==True:
            seventhx= eighthx + 20
        seventhy= eighthy
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy, eighthx, eighthy]
    elif V == True  and mod == ">" and direction == "outie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx + 20
        elif righthanded == True:
            secondx = firstx - 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx - 20
        elif righthanded == True:
            fourthx     = thirdx + 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx - 20
        elif righthanded == True:
            fifthx     = fourthx + 20
        fifthy         = fourthy
        if slant == True and righthanded == False:
            sixthx     =  firstx - 15
        elif slant == True and righthanded == True:
            sixthx   =  firstx + 15
        elif slant == False and righthanded == False:
            sixthx     = fifthx - 20
        elif slant == False and righthanded == True:
            sixthx     = fifthx + 20
        sixthy     = fifthy-25
        if righthanded == False and slant == True:
            eighthx = firstx+10
        elif righthanded==True and slant == True:
            eighthx = firstx-10
        else:
            eighthx = firstx
        eighthy = firsty+20
        if righthanded==False:
            seventhx=eighthx - 20
        elif righthanded==True:
            seventhx= eighthx + 20
        seventhy= eighthy
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy, eighthx, eighthy]
    elif V == True  and mod == ">" and direction == "innie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx - 20
        elif righthanded == True:
            secondx = firstx + 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx + 20
        elif righthanded == True:
            fourthx     = thirdx - 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx + 20
        elif righthanded == True:
            fifthx     = fourthx - 20
        fifthy         = fourthy
        if slant == True:
            sixthx     = fifthx
        elif slant == False and righthanded == False:
            sixthx     = fifthx + 15
        elif slant == False and righthanded == True:
            sixthx     = fifthx - 15
        sixthy     = fifthy-25
        if righthanded == False and slant == True:
            eighthx = firstx+10
        elif righthanded==True and slant == True:
            eighthx = firstx-10
        else:
            eighthx = firstx
        eighthy = firsty+20
        if righthanded==False:
            seventhx=eighthx + 20
        elif righthanded==True:
            seventhx= eighthx - 20
        seventhy= eighthy
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy, eighthx, eighthy]
    elif V == False  and mod == "@" and direction == "outie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx - 20
        elif righthanded == True:
            secondx = firstx + 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx + 20
        elif righthanded == True:
            fourthx     = thirdx - 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx + 20
        elif righthanded == True:
            fifthx     = fourthx - 20
        fifthy         = fourthy
        if slant== True and righthanded == True:
            sixthx     = fourthx+20
        elif slant== True and righthanded == False:
            sixthx     = fourthx-20
        else:
            sixthx     = fourthx
        sixthy     = fifthy-37
        if righthanded == False:
            seventhx   =  firstx + 20
        elif righthanded == True:
            seventhx   =  firstx - 20
        seventhy   =  firsty
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
    elif V == False  and mod == "@" and direction == "innie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx + 20
        elif righthanded == True:
            secondx = firstx - 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx - 20
        elif righthanded == True:
            fourthx     = thirdx + 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx - 20
        elif righthanded == True:
            fifthx     = fourthx + 20
        fifthy         = fourthy
        if slant == True:
            sixthx     = fifthx
        elif slant == False and righthanded == False:
            sixthx     = fifthx + 20
        elif slant == False and righthanded == True:
            sixthx     = fifthx - 20
        sixthy     = fifthy-37
        if righthanded == False:
            seventhx   =  firstx - 20
        elif righthanded == True:
            seventhx   =  firstx + 20
        seventhy   =  firsty
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
    elif V == False  and mod == ">" and direction == "outie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx + 20
        elif righthanded == True:
            secondx = firstx - 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx - 20
        elif righthanded == True:
            fourthx     = thirdx + 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx - 20
        elif righthanded == True:
            fifthx     = fourthx + 20
        fifthy         = fourthy
        if slant == True and righthanded == False:
            sixthx     =  firstx - 20
        elif slant == True and righthanded == True:
            sixthx   =  firstx + 20
        elif slant == False and righthanded == False:
            sixthx     = fifthx + 20
        elif slant == False and righthanded == True:
            sixthx     = fifthx - 20
        sixthy     = fifthy-37
        if righthanded == False:
            seventhx   =  firstx - 20
        elif righthanded == True:
            seventhx   =  firstx + 20
        seventhy   =  firsty
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
    elif V == False  and mod == ">" and direction == "innie" and X == False:
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx = firstx - 20
        elif righthanded == True:
            secondx = firstx + 20
        secondy     = firsty
        if slant == True and righthanded == False:
            thirdx      = secondx + 40
        elif slant == True and righthanded == True:
            thirdx      = secondx - 40
        else:
            thirdx      = secondx
        thirdy          = secondy+ 75
        if righthanded == False:
            fourthx     = thirdx + 20
        elif righthanded == True:
            fourthx     = thirdx - 20
        fourthy         = thirdy
        if righthanded == False:
            fifthx     = fourthx + 20
        elif righthanded == True:
            fifthx     = fourthx - 20
        fifthy         = fourthy
        if slant == True:
            sixthx     = fifthx
        elif slant == False and righthanded == False:
            sixthx     = fifthx + 20
        elif slant == False and righthanded == True:
            sixthx     = fifthx - 20
        sixthy     = fifthy-37
        if righthanded == False:
            seventhx   =  firstx + 20
        elif righthanded == True:
            seventhx   =  firstx - 20
        seventhy   =  firsty
        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]

    elif X == True and mod != "Leucine":
        firstx      = startx
        firsty      = starty
        if righthanded == False:
            secondx     =   startx-30
        elif righthanded== True:
            secondx     =   startx+30
        secondy         =   firsty+37
        if righthanded == False:
            thirdx     =   secondx+30
        elif righthanded== True:
            thirdx     =   secondx-30
        thirdy         =   secondy+ 38
        fourthx         =   thirdx
        fourthy         =   thirdy
        if righthanded == False:
            fifthx      =   fourthx+30
        elif righthanded== True:
            fifthx     =    fourthx-30
        fifthy         =    fourthy-38

        coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy]
    elif mod == "Leucine":
        firstx = startx
        firsty = starty
        secondx= firstx+15
        secondy= firsty+5
        thirdx= secondx-30
        thirdy= secondy+10
        fourthx=thirdx+30
        fourthy=thirdy+10
        fifthx=fourthx-30
        fifthy=fourthy+10
        sixthx=fifthx+30
        sixthy=fifthy+10
        seventhx=sixthx-30
        seventhy=sixthy+10
        eighthx = seventhx+30
        eighthy = seventhy+10
        ninthx  = eighthx-15
        ninthy  = eighthy+5
    #elif mod == "HSA":


        coordinates=[firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy,sixthx,sixthy,seventhx,seventhy,eighthx,eighthy,ninthx,ninthy]
    elif mod == "H":
        firstx, firsty = startx,starty
        secondx, secondy=startx,starty
        thirdx, thirdy=startx,starty
        fourthx, fourthy=startx,starty
        fifthx, fifthy=startx,starty

        coordinates = []
    elif mod == "C":
            firstx      = startx
            firsty      = starty
            if righthanded == False:
                secondx     =   startx-30
            elif righthanded== True:
                secondx     =   startx+30
            secondy         =   firsty+37
            if righthanded == False:
                thirdx     =   secondx+30
            elif righthanded== True:
                thirdx     =   secondx-30
            thirdy         =   secondy+ 38
            fourthx         =   thirdx
            fourthy         =   thirdy
            if righthanded == False:
                fifthx      =   fourthx+30
            elif righthanded== True:
                fifthx     =    fourthx-30
            fifthy         =    fourthy-38
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy]

    if slant == True or mod == "Leucine" or previous_H == True :
        top_bond    = [firstx,firsty+2]
    else:
        top_bond    = [firstx,firsty+20]

    if mod=="Leucine":
        bottom_bond = [ninthx,ninthy]
    else:
        bottom_bond = [fourthx, fourthy-2]

    if slant == True and righthanded == False and mod != "Leucine":
        Labelbond   = [firstx+20, (firsty+thirdy)/2]
    elif slant == True and righthanded == True and mod != "Leucine":
        Labelbond   = [firstx-20, (firsty+thirdy)/2]
    elif mod== "Leucine" and righthanded == False:
        Labelbond   = [firstx-30, (firsty+ninthy)/2]
    elif mod== "Leucine" and righthanded == True:
        Labelbond   = [firstx+30, (firsty+ninthy)/2]
    else:
        Labelbond   = [firstx, (firsty+thirdy)/2]
    return(coordinates, bottom_bond,top_bond,Labelbond)



###############Main programme#######################

root = tk.Tk()

HEIGHT = root.winfo_screenheight()
WIDTH  = root.winfo_screenwidth()
print(HEIGHT,WIDTH)
def browseFiles():
    global textBox
    filename = filedialog.askopenfilename(initialdir = "/",title = "Open File",filetypes = (("Text files","*.txt"),("all files","*.*")))
    entry = Get_input(filename)
    print(entry)
    textBox.delete("1.0","end")
    textBox.insert("1.0",str(entry))
    render_pipeline(lower_canvas)

def save_txt_file():
    to_save = textBox.get("1.0","end")
    f = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    f.write(to_save)
    f.close()

def quit():
    root.destroy()

def donothing():
   filewin = Toplevel(root)
   button = Button(filewin, text="Do nothing button")
   button.pack()
def undo():
    global NOTE_to_redo
    global TYPE_to_redo
    global LENGTH_to_redo
    global MOD_to_redo
    global ANTI_to_redo
    global Polygon_to_redo
    global canvas_polygons
    global canvas_labels
    global temp_label
    global TYPE_labels
    global NOTE_labels
    global MOD_labels
    global ANTI_labels
    global LENGTH_labels
    global deleted_polygons
    global Deletes_to_redo
    label_keyslist = list(canvas_labels.keys())
    if deleted_polygons != {}:
        deleted_polygons_keys = list(deleted_polygons.keys())
        to_replace = max(deleted_polygons_keys)
        domain_coordinates = deleted_polygons.get(to_replace)[0]
        domain_name = deleted_polygons.get(to_replace)[1]
        if "-" not in domain_name:
            if "a" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[0], specificity_colours[1]
            elif "b" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[2], specificity_colours[3]
            elif "c" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[4], specificity_colours[5]
            elif "d" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[6], specificity_colours[7]
            elif "e" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[8], specificity_colours[9]
            elif "f" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[10], specificity_colours[11]
            elif "g" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[12], specificity_colours[13]
            elif "h" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[14], specificity_colours[15]
            elif "X" in str(domain_name):
                heavy_colour, light_colour = specificity_colours[18],specificity_colours[18]
            elif str(domain_name) == "C":
                heavy_colour, light_colour = specificity_colours[19],specificity_colours[19]
            else:
                heavy_colour, light_colour = generic_heavy_colour, generic_light_colour
            if "H" in domain_name:
                domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=heavy_colour, width=2, tags="domain")
            elif "L" in domain_name:
                domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=light_colour, width=2, tags="domain")
            canvas_polygons[domain] = [domain_coordinates, domain_name]
            if Label_lock == True:
                domain_name = re.sub("\.|@|>","",domain_name)
                labelx = get_min_max_coordinates(domain_coordinates)[4]
                labely = get_min_max_coordinates(domain_coordinates)[5]
                label  = lower_canvas.create_text(labelx,labely, text = str(domain_name), tags = "label")
                canvas_labels[label] = [[labelx,labely], domain_name]
        elif domain_name == "-" :####
            domain = lower_canvas.create_line(domain_coordinates, fill=bond_colour, width=2, tags="bonds")
            canvas_polygons[domain] = [domain_coordinates, domain_name]
        elif domain_name == "-H-" :####
            domain = lower_canvas.create_line(domain_coordinates, fill=hinge_colour, width=2, tags="bonds")
            canvas_polygons[domain] = [domain_coordinates, domain_name]
        elif domain_name == "-L-" :####
            domain = lower_canvas.create_line(domain_coordinates, fill=linker_colour, width=2, tags="bonds")
            canvas_polygons[domain] = [domain_coordinates, domain_name]
        elif domain_name == "-disulphide-":####
            domain = lower_canvas.create_line(domain_coordinates, fill=disulphide_colour, width=2, tags=("disulphide","bonds"))
            canvas_polygons[domain] = [domain_coordinates, domain_name]
        deleted_polygons = {}
        Deletes_to_redo[to_replace] = [domain_coordinates, domain_name]
    else:
        keys = list(canvas_polygons.keys())+list(TYPE_labels.keys())+list(NOTE_labels.keys())+list(MOD_labels.keys())+list(ANTI_labels.keys())+list(LENGTH_labels.keys())
        to_delete = max(keys)
        if to_delete in (canvas_polygons.keys()):
            lower_canvas.delete(to_delete)
            Polygon_to_redo[to_delete] = canvas_polygons.get(to_delete)
            domain_coordinates = canvas_polygons.get(to_delete)[0]
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            del canvas_polygons[to_delete]
            for i in range(len(label_keyslist)):
                labelx = canvas_labels.get(label_keyslist[i])[0][0]
                labely = canvas_labels.get(label_keyslist[i])[0][1]
                if x1 <= labelx <= x2 and y1 <= labely <= y2:
                    lower_canvas.delete(label_keyslist[i])
                    del canvas_labels[label_keyslist[i]]
        elif to_delete in list(TYPE_labels.keys()):
            lower_canvas.delete(to_delete)
            TYPE_to_redo[to_delete] = TYPE_labels.get(to_delete)
            del TYPE_labels[to_delete]
        elif to_delete in list(NOTE_labels.keys()):
            lower_canvas.delete(to_delete)
            NOTE_to_redo[to_delete] = NOTE_labels.get(to_delete)
            del NOTE_labels[to_delete]
        elif to_delete in list(MOD_labels.keys()):
            lower_canvas.delete(to_delete)
            MOD_to_redo[to_delete] = MOD_labels.get(to_delete)
            del MOD_labels[to_delete]
        elif to_delete in list(ANTI_labels.keys()):
            lower_canvas.delete(to_delete)
            ANTI_to_redo[to_delete] = ANTI_labels.get(to_delete)
            del ANTI_labels[to_delete]
        elif to_delete in list(LENGTH_labels.keys()):
            lower_canvas.delete(to_delete)
            LENGTH_to_redo[to_delete] = LENGTH_labels.get(to_delete)
            del LENGTH_labels[to_delete]

def redo():
    global Polygon_to_redo
    global NOTE_to_redo
    global TYPE_to_redo
    global LENGTH_to_redo
    global MOD_to_redo
    global ANTI_to_redo
    global canvas_polygons
    global canvas_labels
    global temp_label
    global TYPE_labels
    global NOTE_labels
    global MOD_labels
    global ANTI_labels
    global LENGTH_labels
    global specificity_colours
    global Deletes_to_redo
    global deleted_polygons
    if Deletes_to_redo != {}:
        keys = list(Deletes_to_redo.keys())
        to_delete = max(keys)
        domain_coordinates = Deletes_to_redo.get(to_delete)[0]
        domain_name = Deletes_to_redo.get(to_delete)[1]
        if "-" not in domain_name:
            lower_canvas.delete(to_delete)
            Polygon_to_redo[to_delete] = canvas_polygons.get(to_delete)
            domain_coordinates = canvas_polygons.get(to_delete)[0]
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            del canvas_polygons[to_delete]
            for i in range(len(label_keyslist)):
                labelx = canvas_labels.get(label_keyslist[i])[0][0]
                labely = canvas_labels.get(label_keyslist[i])[0][1]
                if x1 <= labelx <= x2 and y1 <= labely <= y2:
                    lower_canvas.delete(label_keyslist[i])
                    del canvas_labels[label_keyslist[i]]
        elif to_delete in list(TYPE_labels.keys()):
            lower_canvas.delete(to_delete)
            TYPE_to_redo[to_delete] = TYPE_labels.get(to_delete)
            del TYPE_labels[to_delete]
        elif to_delete in list(NOTE_labels.keys()):
            lower_canvas.delete(to_delete)
            NOTE_to_redo[to_delete] = NOTE_labels.get(to_delete)
            del NOTE_labels[to_delete]
        elif to_delete in list(MOD_labels.keys()):
            lower_canvas.delete(to_delete)
            MOD_to_redo[to_delete] = MOD_labels.get(to_delete)
            del MOD_labels[to_delete]
        elif to_delete in list(ANTI_labels.keys()):
            lower_canvas.delete(to_delete)
            ANTI_to_redo[to_delete] = ANTI_labels.get(to_delete)
            del ANTI_labels[to_delete]
        elif to_delete in list(LENGTH_labels.keys()):
            lower_canvas.delete(to_delete)
            LENGTH_to_redo[to_delete] = LENGTH_labels.get(to_delete)
            del LENGTH_labels[to_delete]
        Deletes_to_redo= {}
        deleted_polygons={}
    else:
        keys = list(Polygon_to_redo.keys())+list(NOTE_to_redo.keys())+list(TYPE_to_redo.keys())
        to_redo = min(keys)
        domain_name = Polygon_to_redo.get(to_redo)[1]
        domain_coordinates = Polygon_to_redo.get(to_redo)[0]
        print(domain_name)
        if to_redo in (Polygon_to_redo.keys()):
            if "-" not in domain_name:
                if "a" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[0], specificity_colours[1]
                elif "b" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[2], specificity_colours[3]
                elif "c" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[4], specificity_colours[5]
                elif "d" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[6], specificity_colours[7]
                elif "e" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[8], specificity_colours[9]
                elif "f" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[10], specificity_colours[11]
                elif "g" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[12], specificity_colours[13]
                elif "h" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[14], specificity_colours[15]
                elif "X" in str(domain_name):
                    heavy_colour, light_colour = specificity_colours[18],specificity_colours[18]
                elif str(domain_name) == "C":
                    heavy_colour, light_colour = specificity_colours[19],specificity_colours[19]
                else:
                    heavy_colour, light_colour = generic_heavy_colour, generic_light_colour
                if "H" in domain_name:
                    domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=heavy_colour, width=2, tags="domain")
                elif "L" in domain_name:
                    domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=light_colour, width=2, tags="domain")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
                if Label_lock == True:
                    domain_name = re.sub("\.|@|>","",domain_name)
                    labelx = get_min_max_coordinates(domain_coordinates)[4]
                    labely = get_min_max_coordinates(domain_coordinates)[5]
                    label  = lower_canvas.create_text(labelx,labely, text = str(domain_name), tags = "label")
                    canvas_labels[label] = [[labelx,labely], domain_name]
            elif domain_name == "-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=bond_colour, width=2, tags="bonds")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            elif domain_name == "-H-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=hinge_colour, width=2, tags="bonds")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            elif domain_name == "-L-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=linker_colour, width=2, tags="bonds")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            elif domain_name == "-disulphide-":####
                domain = lower_canvas.create_line(domain_coordinates, fill=disulphide_colour, width=2, tags=("disulphide","bonds"))
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            del Polygon_to_redo[to_redo]
        elif to_redo in list(TYPE_labels.keys()):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "TYPE_labels")
            TYPE_labels[label] = [domain_coordinates, domain_name]
            del TYPE_to_redo[to_redo]
        elif to_redo in list(NOTE_labels.keys()):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "NOTE_labels")
            NOTE_labels[label] = [domain_coordinates, domain_name]
            del NOTE_to_redo[to_redo]
        elif to_redo in list(MOD_labels.keys()):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "MOD_labels")
            MOD_labels[label] = [domain_coordinates, domain_name]
            del MOD_to_redo[to_redo]
        elif to_redo in list(ANTI_labels.keys()):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "ANTI_labels")
            ANTI_labels[label] = [domain_coordinates, domain_name]
            del ANTI_to_redo[to_redo]
        elif to_redo in list(LENGTH_labels.keys()):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "LENGTH_labels")
            LENGTH_labels[label] = [domain_coordinates, domain_name]
            del LENGTH_to_redo[to_redo]







canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH, bg='#E7E0E6')
canvas.pack()

frame = tk.Frame(root, bg = '#80c1ff',bd=5)
frame.place(relwidth=1,relheight=1)

###Input box
textBox = tk.Text(frame, font=40)
textBox.place(relx=0.01,rely = 0.65, relwidth=0.4,relheight=0.15)
###Option box
frame2 = tk.Frame(frame, bg = '#D3D3D3')
frame2.place(relx=0.01, rely = 0.05,relheight = 0.35, relwidth = 0.4)

Domain_Primer = []
Domain_Primer_Lock = ""
domain_type = ""
domain_mod  = ""
extra_mods  = ""
domain_charge=""
domain_direction = "constant"
a_button= tk.Button(frame2,text="a",bg = "grey", command =lambda: domain_type_button("a"))
a_button.place(relx = 0.61, rely = 0.21, relheight = 0.1, relwidth=0.05)
b_button= tk.Button(frame2,text="b",bg = "grey", command =lambda: domain_type_button("b"))
b_button.place(relx = 0.66, rely = 0.21, relheight = 0.1, relwidth=0.05)
c_button= tk.Button(frame2,text="c",bg = "grey", command =lambda: domain_type_button("c"))
c_button.place(relx = 0.71, rely = 0.21, relheight = 0.1, relwidth=0.05)
d_button= tk.Button(frame2,text="d",bg = "grey", command =lambda: domain_type_button("d"))
d_button.place(relx = 0.76, rely = 0.21, relheight = 0.1, relwidth=0.05)
e_button= tk.Button(frame2,text="e",bg = "grey", command =lambda: domain_type_button("e"))
e_button.place(relx = 0.61, rely = 0.31, relheight = 0.1, relwidth=0.05)
f_button= tk.Button(frame2,text="f",bg = "grey", command =lambda: domain_type_button("f"))
f_button.place(relx = 0.66, rely = 0.31, relheight = 0.1, relwidth=0.05)
g_button= tk.Button(frame2,text="g",bg = "grey", command =lambda: domain_type_button("g"))
g_button.place(relx = 0.71, rely = 0.31, relheight = 0.1, relwidth=0.05)
h_button= tk.Button(frame2,text="h",bg = "grey", command =lambda: domain_type_button("h"))
h_button.place(relx = 0.76, rely = 0.31, relheight = 0.1, relwidth=0.05)
Mod_button = tk.Button(frame2,text="*",bg = "grey", command =lambda: extra_mod_button("*"))
Mod_button.place(relx = 0.61, rely = 0.41, relheight = 0.1, relwidth=0.05)
Gly_button = tk.Button(frame2,text="!",bg = "grey", command =lambda: extra_mod_button("!"))
Gly_button.place(relx = 0.66, rely = 0.41, relheight = 0.1, relwidth=0.05)
KIH_knob= tk.Button(frame2,text=">",bg = "grey", command =lambda: domain_mod_button(">"))
KIH_knob.place(relx = 0.71, rely = 0.41, relheight = 0.1, relwidth=0.05)
KIH_hole= tk.Button(frame2,text="@",bg = "grey", command =lambda: domain_mod_button("@"))
KIH_hole.place(relx = 0.76, rely = 0.41, relheight = 0.1, relwidth=0.05)
Positive_charge= tk.Button(frame2,text="+",bg = "grey", command =lambda: domain_charge_button("+"))
Positive_charge.place(relx = 0.61, rely = 0.51, relheight = 0.1, relwidth=0.05)
Negative_charge= tk.Button(frame2,text="-",bg = "grey", command =lambda: domain_charge_button("_"))
Negative_charge.place(relx = 0.66, rely = 0.51, relheight = 0.1, relwidth=0.05)
LengthLabelButton = tk.Button(frame2,text="LENGTH", bg="grey", command=lambda: SelectCommentTypeButton("LENGTH"))
LengthLabelButton.place(relx = 0.71,rely = 0.51, relheight = 0.1, relwidth= 0.1)
NoteLabelButton = tk.Button(frame2,text="NOTE", bg="grey", command=lambda: SelectCommentTypeButton("NOTE"))
NoteLabelButton.place(relx = 0.61,rely = 0.61, relheight = 0.1, relwidth= 0.1)
TypeLabelButton = tk.Button(frame2,text="TYPE", bg="grey", command=lambda: SelectCommentTypeButton("TYPE"))
TypeLabelButton.place(relx = 0.71,rely = 0.61, relheight = 0.1, relwidth= 0.1)
AntiLabelButton = tk.Button(frame2,text="ANTI", bg="grey", command=lambda: SelectCommentTypeButton("ANTI"))
AntiLabelButton.place(relx = 0.61,rely = 0.71, relheight = 0.1, relwidth= 0.1)
ModLabelButton = tk.Button(frame2,text="MOD", bg="grey", command=lambda: SelectCommentTypeButton("MOD"))
ModLabelButton.place(relx = 0.71,rely = 0.71, relheight = 0.1, relwidth= 0.1)
nanobody_button = tk.Button(frame2,text="nanobody",bg = "grey", command =lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"Single_Fv_Chain",False,domain_mod,"","",str("nanoVH"+domain_mod+"."+domain_type),False,True))
nanobody_button.place(relx = 0.61, rely = 0.01, relheight = 0.2, relwidth=0.2)
###Insert bonds buttons ###
##Col1
InsertVHDomainButton= tk.Button(frame2,text="VH",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"innie",False,domain_mod,"","",str("VH"+domain_mod+"."+domain_type),False,True))
InsertVHDomainButton.place(relx = 0.21, rely = 0.01, relheight = 0.2, relwidth=0.2)
InsertCH1DomainButton= tk.Button(frame2,text="CH1",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH1"+domain_mod),False,True))
InsertCH1DomainButton.place(relx = 0.21, rely = 0.21, relheight = 0.2, relwidth=0.2)
InsertCH2DomainButton= tk.Button(frame2,text="CH2",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH2"+domain_mod),False,True))
InsertCH2DomainButton.place(relx = 0.21, rely = 0.41, relheight = 0.2, relwidth=0.2)
InsertCH3DomainButton= tk.Button(frame2,text="CH3",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH3"+domain_mod),False,True))
InsertCH3DomainButton.place(relx = 0.21, rely = 0.61, relheight = 0.2, relwidth=0.2)
##Col2
InsertVLDomainButton= tk.Button(frame2,text="VL",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"outie",False,domain_mod,"","",str("VL"+domain_mod+"."+domain_type),True,False))
InsertVLDomainButton.place(relx = 0.41, rely = 0.01, relheight = 0.2, relwidth=0.2)
InsertCLDomainButton= tk.Button(frame2,text="CL",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CL"+domain_mod),True,False))
InsertCLDomainButton.place(relx = 0.41, rely = 0.21, relheight = 0.2, relwidth=0.2)
InsertCH4DomainButton= tk.Button(frame2,text="CH4",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH4"+domain_mod),False,True))
InsertCH4DomainButton.place(relx = 0.41, rely = 0.41, relheight = 0.2, relwidth=0.2)
InsertXDomainButton= tk.Button(frame2,text="Other",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,True,domain_mod,"","",str("X"+domain_mod),True,False))
InsertXDomainButton.place(relx = 0.41, rely = 0.61, relheight = 0.2, relwidth=0.1)
InsertCDomainButton= tk.Button(frame2,text="Chem",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,True,domain_mod,"","",str("C"+domain_mod),True,False))
InsertCDomainButton.place(relx = 0.51, rely = 0.61, relheight = 0.2, relwidth=0.1)

###Drag and pull bonds###
Bond_lock = ""
InsertBondButton= tk.Button(frame2,text="Connect",bg = "grey", command=lambda: bond_drag_button(lower_canvas,"-","bond"))
InsertBondButton.place(relx = 0.01, rely = 0.01, relheight = 0.2, relwidth=0.2)
InsertDBondButton= tk.Button(frame2,text="Disulphide",bg = "grey", command=lambda: disulphide_bond_button(lower_canvas,"-","disulphide"))
InsertDBondButton.place(relx = 0.01, rely = 0.21, relheight = 0.2, relwidth=0.2)
InsertLHingeButton= tk.Button(frame2,text="Hinge",bg = "grey", command=lambda: Hinge_bond_button(lower_canvas,"-H-","hinge"))
InsertLHingeButton.place(relx = 0.01, rely = 0.41, relheight = 0.2, relwidth=0.2)
InsertLinkerButton= tk.Button(frame2,text="Linker",bg = "grey", command=lambda: Linker_bond_button(lower_canvas,"-L-","linker"))
InsertLinkerButton.place(relx = 0.01, rely = 0.61, relheight = 0.2, relwidth=0.2)

###Delete/clear###
Delete_lock = False
Label_lock = True
CustomLabelLock = ""
InsertDelAllButton = tk.Button(frame2,text="Clear All",bg="grey", command=lambda: delete_all_button(lower_canvas))
InsertDelAllButton.place(relx = 0.81,rely = 0.01, relheight = 0.2, relwidth= 0.18)
InsertDelClickButton = tk.Button(frame2,text="Delete",bg="grey", command=lambda: delete_button(lower_canvas))
InsertDelClickButton.place(relx = 0.81,rely = 0.21, relheight = 0.2, relwidth= 0.18)
Labels_buttons = tk.Button(frame2,text="Labels", bg="grey", command=lambda: labels_button(lower_canvas))
Labels_buttons.place(relx = 0.81,rely = 0.41, relheight = 0.2, relwidth= 0.18)
CustomLabelButton = tk.Button(frame2,text="Comment", bg="grey", command=lambda: CommentLabelButton_function(lower_canvas))
CustomLabelButton.place(relx = 0.81,rely = 0.61, relheight = 0.2, relwidth= 0.18)
CustomLabelEntry = tk.Text(frame2, font=40)
CustomLabelEntry.place(relx=0.01,rely = 0.83, relwidth=0.98,relheight=0.15)
##Status bar###
status_label = tk.Label(root, text='test', bd=1)
status_label.place(rely = 0.98, relheight = 0.02, relwidth = 1)

##Big button
button = tk.Button(frame, text = "Get Structure", bg = "grey", font=40, command=lambda: render_pipeline(lower_canvas))
button.place(relx=0.05,rely=0.825,relheight=0.1, relwidth=0.15)
button.bind("<Enter>", button_hover)
button.bind("<Leave>", button_hover_leave)
button = tk.Button(frame, text = "Get Sequence", bg = "grey", font=40, command=lambda: sequence_pipeline(lower_canvas))
button.place(relx=0.2,rely=0.825,relheight=0.1, relwidth=0.15)
button.bind("<Enter>", button_hover)
button.bind("<Leave>", button_hover_leave)
button = tk.Button(frame, text = "Tidy", bg = "grey", font=40, command=lambda: sequence_render_pipeline(lower_canvas))
button.place(relx=0.05,rely=0.925,relheight=0.05, relwidth=0.3)
button.bind("<Enter>", button_hover)
button.bind("<Leave>", button_hover_leave)


###Library
#frame1 = tk.Frame(frame, bg = '#D3D3D3')
#frame1.place(relx=0.01, rely = 0.425,relheight = 0.2, relwidth = 0.4)
Library= tk.Listbox(frame, selectbackground='#D3D3D3', height=20)
antibodyformats = {
"IgG":"VH.a(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3(5:12) | VL.a(6:1)-CL(7:2){1} | VH.a(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3(12:5) | VL.a(13:8)-CL(14:9){1}",
"Promiscuous IgG": "VH.ab(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3(5:12) | VL.ab(6:1)-CL(7:2){1} | VH.cd(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3(12:5) | VL.cd(13:8)-CL(14:9){1}",
"Knobs in Holes":"VH.a(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3@(5:12) | VL.a(6:1)-CL(7:2){1} | VH.b(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3>(12:5) | VL.b(13:8)-CL(14:9){1}",
"Orthagonal Fab":"VH+.a(1:6)-CH1>(2:7){1}-H(3:10){2}-CH2(4:11)-CH3@(5:12)|VL_.a(6:1)-CL@(7:2){1}|VH.b(8:13)-CH1>(9:14){1}-H(10:3){2}-CH2(11:4)-CH3>(12:5)|VL.b(13:8)-CL@(14:9){1}",
"IgG-H-scFV":"VH.a(1:8)-CH1(2:9){1}-H(3:12){2}-CH2(4:13)-CH3(5:14)-L-VH.b(6:7)-L-VL.b(7:6)|VL.a(8:1)-CL(9:2){1}|VH.a(10:17)-CH1(11:18){1}-H(12:3){2}-CH2(13:4)-CH3(14:5)-L-VH.b(15:16)-L-VL.b(16:15)|VL.a(17:10)-CL(18:11){1}",
"IgG-L-scFV":"VH.a(1:6) -CH1(2:7){1} -H(3:12){2} -CH2(4:13) -CH3(5:14) |VL.a(6:1) -CL(7:2){1} -L -VH.b(8:9) -L -VL.b(9:8) |VH.a(10:15) -CH1(11:16){1} -H(12:3){2} -CH2(13:4) -CH3(14:5) | VL.a(15:10) -CL(16:11){1} -L -VH.b(17:18) -L -VL.b(18:17)",
"scFV-H-scFV": "VL.a(1:2)-L-VH.a(2:1)-H(3:8){2}-VH.b(4:5)-L-VL.b(5:4)|VL.c(6:7)-L-VH.c(7:6)-H(8:3){2}-VH.d(9:10)-L-VL.d(10:9)",
"F(ab)2":"VH.a(1:4) - CH1(2:5){1} -H(3:8){2} | VL.a(4:1) - CL(5:2){1} | VH.b(6:9) - CH1(7:10){1} - H(8:3){2} | VL.b(9:6) - CL(10:7){1}",
"scFV":"VH.a(1:2)-L-VL.a(2:1)",
"scFV4":"VL.a(4:5)-L-VH.a(5:4)-CH1(6:3){1}-H(7:16){2}-CH2(8:17)-CH3(9:18)|  VH.a(1:2)-L-VL.a(2:1)-CL(3:6){1} | VL.b(13:14)-L-VH.b(14:13)-CH1(15:12){1}-H(16:7){2}-CH2(17:8)-CH3(18:9) | VH.b(10:11)-L-VL.b(11:10)-CL(12:15){1}",
"sdFV4":"VH.a(1) -CH1(2:7){1} -H(3:10){2}-CH2(4:11) -CH3@(5:12) | VH.a(6)  -CL(7:2){1} | VH.b(8) -CH1(9:14){1}-H(10:3){2}-CH2(11:4) -CH3>(12:5) | VH.b(13) -CL(14:9){1}",
"Nanobody":"VH.a(1)",
"BiTE":"VH.a(1:2) -L -VL.a(2:1) -L -VH.b(3:4) -L -VL.b(4:3)",
"HSAbody":"VL.a(1:2)-L-VH.a(2:1)-L-X(3)[TYPE:FUSION][NOTE:human serum albumin]-L-VH.b(4:5)-L-VL.b(5:4)",
"Cov-X-body":"X(1)[TYPE: pharmacophore peptide heterodimer]-VH.a(2:7)- CH1(3:8){1}-H(4:12){2}-CH2(5:11)-CH3(6:12) | VL.a(7:2)-CL(8:3){1} | X(9)[TYPE: pharmacophore peptide heterodimer]-VH.b(10:15)-CH1(11:16){1}-H(12:4){2}- CH2(13:5)-CH3 (14:6) | VL.b(15:10)-CL(16:11){1}",
"Diabody":"VH.a(1:3) -L -VH.b(2:4) | VL.a(3:1) -L -VL.b(4:2)",
"Miniantibody":"VH.a(1:2)-L-VL.a(2:1)-H*(3:7){1}[MOD: engineered disulphide bond]-X(4)[TYPE:LEUCINE] | VH.b(5:6)-L-VL.b(6:5)-H*(7:3){1}[MOD: engineered disulphide bond]-X(8)[TYPE:LEUCINE]",
"scDiabody":"VH.a(1:4) -L -VL.a(2:3) -L -VH.b(3:2) -L -VL.b(4:1)",
"scDiabody-CH3":"VL.a(1:4){1} -L -VL.a(2:3) -L -VH.a(3:2) -L -VH.a(4:1){1}-H(5:12){2}-CH2(6:13)-CH3(7:14)| VL.b(8:11){1}-L-VL.b(9:10)-L-VH.b(10:9)-L-VH.b(11:8){1}-H(12:5){2}-CH2(13:6)-CH3(14:7)",
"scDiabody-Fc":"VH.b(1:2) -L -VL.b(2:1) -L -VH.a(3:4) -L -VL.a(4:3) -H(5:12){2} -CH2(6:13) -CH3(7:14) | VH.b(8:9) -L -VL.b(9:8) -L -VH.a(10:11) -L -VL.a(11:10) -H(12:5){2} -CH2(13:6) -CH3(14:7)",
"DART":"VL.a(1:5) -L -VH.b(2:4) -H*(3:6){3}[MOD: engineered disulphide bond] | VL.b(4:2) -L -VH.a(5:1) -H*(6:3){3}[MOD: engineered disulphide bond]",
"Tandem A and B": "VH.a(1:5) -L -VL.b(2:6) -L -VH.b(3:7) -L -VL.a(4:8) | VL.a(5:1) -L -VH.b(6:2) -L -VL.b(7:3) -L -VH.a(8:4)",
"Intrabody":"VL.a(1:2) -L -VH.a(2:1)-H(3:10){2}- CH1(4:11) - CH2(5:12) -L - VH.b(6:7) -L -VL.b(7:6)| VL.a(8:9) -L -VH.a(9:8) -H(10:3){2} -CH1(11:4) - CH2(12:5) -L -VL.b(13:14) -L -VH.b(14:13)",
"Fv-Fc":"VH.a(1:5){1}-H(2:7){2}-CH1(3:8)-CH2(4:9)| VL.a(5:1){1}|VH.b(6:10){1}-H(7:2){2}-CH1(8:3)-CH2(9:4)|VL.b(10:6){1}",
"Triplebody":"VH.a(1:5)-CH1(2:6){1} -L- VL.b(3:4) -L -VH.b(4:3) | VL.a(5:1) -CL(6:2){1} -L - VL.c(7:8) -L -VH.c(8:7)",
"scTriplebody":"VH.a(1:5)-CH1(2:6){2}-L-VL.b(3:4){1}-L-VH.b(4:3){1}-L-VL.a(5:1)-CL(6:2){2}-L-VL.c(7:8){1}-L-VH.c(8:7){1}",
"TriBiMinibody":"VH.a(1:2) -L -VL.a(2:1) -H(3:9){2} -CH3@(4:10){2} -L -VH.b(5:6) -L - VL.b(6:5) | VH.c(7:8) -L -VL.c(8:7) -H(9:3){2}-CH3>(10:4)",
"LUZ-Y":"VL.a(1:3)-CL(2:4){1}-L -VH.a(3:1)-CH1(4:2){1}-H(5:13){2}-CH2(6:14) -CH3(7:15) -X(8)[TYPE: LEUCINE] | VL.b(9:11)-CL(10:12){1} -L -VH.b(11:9) -CH1(12:10){1}-H(13:5){2}-CH2(14:6)-CH3(15:7)-X(16)[TYPE: LEUCINE]",
"Dock and Lock":"VH.a(1:4)-CH1(2:5){1}-L-X(3)[TYPE:FUSION]|VL.a(4:1)-CL(5:2){1}|VH.b(6:8)-CH1(7:9){1}-L-X(3)|VL.b(8:6)-CL(9:7){1}|VH.c(10:12)-CH1(11:13){1}-L-X(3)|VL.c(12:10)-CL(13:11){1}|VH.d(14:16)-CH1(15:17){1}-L-X(3)|VL.d(16:14)-CL(17:15){1}",
"scFV-IgG-scFV-scFV": "VL.b(1:2)-L-VH.b(2:1)-L-VH.a(3:12)-CH1(4:13){1}-H(5:18){2}-CH2(6:19)-CH3(7:20)-L-VH.b(8:9)-L-VL.b(9:8)-L-VH.c(10:11)-L-VL.c(11:10)|VL.a(12:3)-CL(13:4){1}|VL.b(14:15)-L-VH.b(15:14)-L-VH.a(16:25)-CH1(17:26){1}-H(18:5){2}-CH2(19:6)-CH3(20:7)-L-VH.b(21:22)-L-VL.b(22:21)-L-VH.c(23:24)-L-VL.c(24:23)|VL.a(25:16)-CL(26:17){1}",
"scFV-scFV-Fc":"VH.a(1:2)-L-VL.a(2:1)-L-VH.b(3:4)-L-VL.b(4:3)-CH2(5:7)-CH3(6:8)-L-CH2(7:5)-CH3(8:6)",
"Trimeric Fusion Protein":"VH.a(1:6)-CH1(2:7){1}-H(3:11){2}-CH2(4:12)-CH3(5:13)|VL.a(6:1)-CL(7:2){1}|X(8:9,14)[NOTE:FUSION]-X(9:8,14)[NOTE:FUSION]-CH1(10:15){1}-H(11:3){2}-CH2(12:4)-CH3(13:5)|X(14:8,9)[NOTE:FUSION]-CL(15:10){1}",
"IgG-IgG":"VH.a(1:12)-CH1(2:13){1}-H(3:8){2}-CH2(4:9)-CH3(5:10)|VL.a(12:1)-CL(13:2){1}|VH.a(6:14)-CH1(7:15){1}-H(8:3){2}-CH2(9:4)-CH3(10:5)-L-C(11)[MOD: orthophenylenedimaleimide fusion]|VL.a(14:6)-CL(15:7){1}|VH.b(16:26)-CH1(17:27){1}-H(18:23){2}-CH2(19:24)-CH3(20:25)-L-C(11)[MOD: orthophenylenedimaleimide fusion]|VH.b(21:28)-CH1(22:29){1}-H(23:18){2}-CH2(24:19)-CH3(25:20)|VL.b(26:16)-CL(27:17){1}|VL.b(28:21)-CL(29:22){1}"
}

formats_keyslist= list(antibodyformats.keys())
for i in range(len(formats_keyslist)):
    Library.insert("end",formats_keyslist[i])
    if i %2==0:
        Library.itemconfig(i, bg='#D3D3D3')


Library.place(relx=0.01, rely = 0.425,relheight = 0.2, relwidth = 0.4)
Library.bind('<<ListboxSelect>>', items_selected)


#Library.place(frame1, relheight=0.05, relwidth=0.3)
####Colours
a_heavy_colour, a_light_colour = '#007ECB', '#73CAFF'
b_heavy_colour, b_light_colour = '#FF43EE', '#F9D3F5'
c_heavy_colour, c_light_colour = '#0BD05A', '#B9FAD3'
d_heavy_colour, d_light_colour = '#D9DE4A', '#E2F562'
e_heavy_colour, e_light_colour = '#FF0000', '#FF8585'
f_heavy_colour, f_light_colour = '#FE6D03', '#FFC296'
g_heavy_colour, g_light_colour = '#5C54FF', '#AFABFF'
h_heavy_colour, h_light_colour = '#0FFBFF', '#CCFEFF'
generic_heavy_colour, generic_light_colour = '#8E8E8E','#C6C4C4'
disulphide_colour = "red"
bond_colour = "black"
hinge_colour = "dark green"
linker_colour = "purple"
X_colour = "#68C1C1"
C_colour = "#ABA600"

specificity_colours = [a_heavy_colour, a_light_colour,b_heavy_colour, b_light_colour,c_heavy_colour, c_light_colour,d_heavy_colour, d_light_colour,e_heavy_colour, e_light_colour,f_heavy_colour, f_light_colour,g_heavy_colour, g_light_colour,h_heavy_colour, h_light_colour,generic_heavy_colour,generic_light_colour,X_colour,C_colour]
bond_colours = [disulphide_colour,bond_colour,hinge_colour,linker_colour]

#Coloursettings = tk.Listbox(frame, selectbackground='#D3D3D3', height=20)
def open_settings():
    top = tk.Toplevel()
    top.title = "Settings"
    top.geometry('500x400')
    tabControl = ttk.Notebook(top)

    tab1 = ttk.Frame(tabControl)
    tab2 = ttk.Frame(tabControl)

    tabControl.add(tab1, text ='Pairing Sensitivity')
    tabControl.add(tab2, text ='Colour Changer')
    tabControl.pack(expand = 1, fill ="both")

    ttk.Label(tab1,text ="Pairing").grid(column = 0,
                                   row = 0,
                                   padx = 30,
                                   pady = 30)

    ttk.Label(tab2,text ="Domains").place(relx= 0.1, rely = 0.1)
    ttk.Label(tab2,text = "Colour changer")
    global specificity_colours
    global bond_colours
    selected_domain = ""
    selected_colour = ""
    colourindex = 0
    def recolorise(string1,string2,new_colour):
        global canvas_polygons
        polygons_keyslist = list(canvas_polygons.keys())
        for i in range(len(polygons_keyslist)):
            if string1 in canvas_polygons.get(polygons_keyslist[i])[1] and string2 in canvas_polygons.get(polygons_keyslist[i])[1]:
                lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
    def browse_colour():
        global selected_colour
        global colourindex
        global selected_domain
        global canvas_polygons
        new_colour = colorchooser.askcolor()[1]
        selected_colour = new_colour
        colours[coloursettings_keyslist[colourindex]] = new_colour
        polygons_keyslist = list(canvas_polygons.keys())

        if   selected_domain == "VHa":
            specificity_colours[0] = new_colour
            recolorise("VH","a",new_colour)
        elif selected_domain == "VLa":
            specificity_colours[1] = new_colour
            recolorise("VL","a",new_colour)
        elif selected_domain == "VHb":
            specificity_colours[2] = new_colour
            recolorise("VH","b",new_colour)
        elif selected_domain == "VLb":
            specificity_colours[3] = new_colour
            recolorise("VL","b",new_colour)
        elif selected_domain == "VHc":
            specificity_colours[4] = new_colour
            recolorise("VH","c",new_colour)
        elif selected_domain == "VLc":
            specificity_colours[5] = new_colour
            recolorise("VL","c",new_colour)
        elif selected_domain == "VHd":
            specificity_colours[6] = new_colour
            recolorise("VH","d",new_colour)
        elif selected_domain == "VLd":
            specificity_colours[7] = new_colour
            recolorise("VL","d",new_colour)
        elif selected_domain == "VHe":
            specificity_colours[8] = new_colour
            recolorise("VH","e",new_colour)
        elif selected_domain == "VLe":
            specificity_colours[9] = new_colour
            recolorise("VL","e",new_colour)
        elif selected_domain == "VHf":
            specificity_colours[10] = new_colour
            recolorise("VH","f",new_colour)
        elif selected_domain == "VLf":
            specificity_colours[11] = new_colour
            recolorise("VL","f",new_colour)
        elif selected_domain == "VHg":
            specificity_colours[12] = new_colour
            recolorise("VH","g",new_colour)
        elif selected_domain == "VLg":
            specificity_colours[13] = new_colour
            recolorise("VL","g",new_colour)
        elif selected_domain == "VHh":
            specificity_colours[14] = new_colour
            recolorise("VH","h",new_colour)
        elif selected_domain == "VLh":
            specificity_colours[15] = new_colour
            recolorise("VL","h",new_colour)
        elif selected_domain == "X":
            specificity_colours[18] = new_colour
        elif selected_domain == "C":
            specificity_colours[19] = new_colour
            for i in range(len(polygons_keyslist)):
                if "X" in canvas_polygons.get(polygons_keyslist[i])[1]:
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        elif selected_domain == "-":
            bond_colours[1] = new_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "-":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        elif selected_domain == "-L-":
            bond_colours[3] = new_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "-L-":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        elif selected_domain == "-H-":
            bond_colours[2] = new_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "-H-":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        elif selected_domain == "Disulphide":
            bond_colours[0] = new_colour
            disulphides_list = lower_canvas.find_withtag("disulphide")
            disulphides_dict = {}
            for i in range(len(disulphides_list)):
                for j in range(len(polygons_keyslist)):
                    if disulphides_list[i] == polygons_keyslist[j]:
                        disulphides_dict[j] = canvas_polygons.get(polygons_keyslist[j])
            disulphides_keyslist = list(disulphides_dict.keys())
            for i in range(len(disulphides_keyslist)):
                lower_canvas.itemconfig(disulphides_keyslist[i], fill = new_colour)
        colour_checker.config(bg=new_colour)

    def revertcolor():
        global specificity_colours
        global bond_colours
        global selected_domain
        global canvas_polygons
        polygons_keyslist = list(canvas_polygons.keys())
        if   selected_domain == "VHa":
            specificity_colours[0] = a_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = a_heavy_colour
            recolorise("VH","a",a_heavy_colour)
            colour_checker.config(bg=a_heavy_colour)
        elif selected_domain == "VLa":
            specificity_colours[1] = a_light_colour
            colours[coloursettings_keyslist[colourindex]] = a_light_colour
            recolorise("VL","a",a_light_colour)
            colour_checker.config(bg=a_light_colour)
        elif selected_domain == "VHb":
            specificity_colours[2] = b_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = b_heavy_colour
            recolorise("VH","b",b_heavy_colour)
            colour_checker.config(bg=b_heavy_colour)
        elif selected_domain == "VLb":
            specificity_colours[3] = b_light_colour
            colours[coloursettings_keyslist[colourindex]] = b_light_colour
            recolorise("VL","b",b_light_colour)
            colour_checker.config(bg=b_light_colour)
        elif selected_domain == "VHc":
            specificity_colours[4] = c_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = c_heavy_colour
            recolorise("VH","c",c_heavy_colour)
            colour_checker.config(bg=c_heavy_colour)
        elif selected_domain == "VLc":
            specificity_colours[5] = c_light_colour
            colours[coloursettings_keyslist[colourindex]] = c_light_colour
            recolorise("VL","c",c_light_colour)
            colour_checker.config(bg=c_light_colour)
        elif selected_domain == "VHd":
            specificity_colours[6] = d_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = d_heavy_colour
            recolorise("VH","d",d_heavy_colour)
            colour_checker.config(bg=d_heavy_colour)
        elif selected_domain == "VLd":
            specificity_colours[7] = d_light_colour
            colours[coloursettings_keyslist[colourindex]] = d_light_colour
            recolorise("VL","d",d_light_colour)
            colour_checker.config(bg=d_light_colour)
        elif selected_domain == "VHe":
            specificity_colours[8] = e_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = e_heavy_colour
            recolorise("VH","e",e_heavy_colour)
            colour_checker.config(bg=e_heavy_colour)
        elif selected_domain == "VLe":
            specificity_colours[9] = e_light_colour
            colours[coloursettings_keyslist[colourindex]] = e_light_colour
            recolorise("VL","e",e_light_colour)
            colour_checker.config(bg=e_light_colour)
        elif selected_domain == "VHf":
            specificity_colours[10] = f_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = f_heavy_colour
            recolorise("VH","f",f_heavy_colour)
            colour_checker.config(bg=f_heavy_colour)
        elif selected_domain == "VLf":
            specificity_colours[11] = f_light_colour
            colours[coloursettings_keyslist[colourindex]] = f_light_colour
            recolorise("VL","f",f_light_colour)
            colour_checker.config(bg=f_light_colour)
        elif selected_domain == "VHg":
            specificity_colours[12] = g_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = g_heavy_colour
            recolorise("VH","g",g_heavy_colour)
            colour_checker.config(bg=g_heavy_colour)
        elif selected_domain == "VLg":
            specificity_colours[13] = g_light_colour
            colours[coloursettings_keyslist[colourindex]] = g_light_colour
            recolorise("VL","g",g_light_colour)
            colour_checker.config(bg=g_light_colour)
        elif selected_domain == "VHh":
            specificity_colours[14] = h_heavy_colour
            colours[coloursettings_keyslist[colourindex]] = h_heavy_colour
            recolorise("VH","h",h_heavy_colour)
            colour_checker.config(bg=h_heavy_colour)
        elif selected_domain == "VLh":
            specificity_colours[15] = h_light_colour
            colours[coloursettings_keyslist[colourindex]] = h_light_colour
            recolorise("VL","h",h_light_colour)
            colour_checker.config(bg=h_light_colour)
        elif selected_domain == "X":
            specificity_colours[18] = X_colour
            for i in range(len(polygons_keyslist)):
                if "X" in canvas_polygons.get(polygons_keyslist[i])[1]:
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = X_colour)
            colours[coloursettings_keyslist[colourindex]] = X_colour
            colour_checker.config(bg=X_colour)
        elif selected_domain == "C":
            specificity_colours[19] = C_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "C":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = C_colour)
            colours[coloursettings_keyslist[colourindex]] = C_colour
            colour_checker.config(bg=C_colour)
        elif selected_domain == "-":
            bond_colours[1] = bond_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "-":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = bond_colour)
            colours[coloursettings_keyslist[colourindex]] = bond_colour
            colour_checker.config(bg=bond_colour)
        elif selected_domain == "-L-":
            bond_colours[3] = linker_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "-L-":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = linker_colour)
            colours[coloursettings_keyslist[colourindex]] = linker_colour
            colour_checker.config(bg=linker_colour)
        elif selected_domain == "-H-":
            bond_colours[2] = hinge_colour
            for i in range(len(polygons_keyslist)):
                if canvas_polygons.get(polygons_keyslist[i])[1] == "-H-":
                    lower_canvas.itemconfig(polygons_keyslist[i], fill = hinge_colour)
            colours[coloursettings_keyslist[colourindex]] = hinge_colour
            colour_checker.config(bg=hinge_colour)
        elif selected_domain == "Disulphide":
            bond_colours[0] = disulphide_colour
            disulphides_list = lower_canvas.find_withtag("disulphide")
            disulphides_dict = {}
            for i in range(len(disulphides_list)):
                for j in range(len(polygons_keyslist)):
                    if disulphides_list[i] == polygons_keyslist[j]:
                        disulphides_dict[j] = canvas_polygons.get(polygons_keyslist[j])
            disulphides_keyslist = list(disulphides_dict.keys())
            for i in range(len(disulphides_keyslist)):
                lower_canvas.itemconfig(disulphides_keyslist[i], fill = disulphide_colour)
            colours[coloursettings_keyslist[colourindex]] = disulphide_colour
            colour_checker.config(bg=disulphide_colour)




    def revertallcolours():
        global specificity_colours
        global bond_colours
        global selected_domain
        global selected_colour
        specificity_colours[0] = a_heavy_colour
        specificity_colours[1] = a_light_colour
        specificity_colours[2] = b_heavy_colour
        specificity_colours[3] = b_light_colour
        specificity_colours[4] = c_heavy_colour
        specificity_colours[5] = c_light_colour
        specificity_colours[6] = d_heavy_colour
        specificity_colours[7] = d_light_colour
        specificity_colours[8] = e_heavy_colour
        specificity_colours[9] = e_light_colour
        specificity_colours[10] = f_heavy_colour
        specificity_colours[11] = f_light_colour
        specificity_colours[12] = g_heavy_colour
        specificity_colours[13] = g_light_colour
        specificity_colours[14] = h_heavy_colour
        specificity_colours[15] = h_light_colour
        specificity_colours[18] = X_colour
        specificity_colours[19] = C_colour
        bond_colours[1] = bond_colour
        bond_colours[3] = linker_colour
        bond_colours[2] = hinge_colour
        bond_colours[0] = disulphide_colour
        recolorise("VH","a",a_heavy_colour)
        recolorise("VH","b",b_heavy_colour)
        recolorise("VH","c",c_heavy_colour)
        recolorise("VH","d",d_heavy_colour)
        recolorise("VH","e",e_heavy_colour)
        recolorise("VH","f",f_heavy_colour)
        recolorise("VH","g",g_heavy_colour)
        recolorise("VH","h",h_heavy_colour)
        recolorise("VL","a",a_light_colour)
        recolorise("VL","b",b_light_colour)
        recolorise("VL","c",c_light_colour)
        recolorise("VL","d",d_light_colour)
        recolorise("VL","e",e_light_colour)
        recolorise("VL","f",f_light_colour)
        recolorise("VL","g",g_light_colour)
        recolorise("VL","h",h_light_colour)
        colours["VHa"] = a_heavy_colour
        colours["VHb"] = b_heavy_colour
        colours["VHc"] = c_heavy_colour
        colours["VHd"] = d_heavy_colour
        colours["VHe"] = e_heavy_colour
        colours["VHf"] = f_heavy_colour
        colours["VHg"] = g_heavy_colour
        colours["VHh"] = h_heavy_colour
        colours["VLa"] = a_light_colour
        colours["VLb"] = b_light_colour
        colours["VLc"] = c_light_colour
        colours["VLd"] = d_light_colour
        colours["VLe"] = e_light_colour
        colours["VLf"] = f_light_colour
        colours["VLg"] = g_light_colour
        colours["VLh"] = h_light_colour
        colours["X"]   = X_colour
        colours["C"]   = C_colour
        colours["-"]   = bond_colour
        colours["-L-"] = linker_colour
        colours["-H-"] = hinge_colour
        colours["Disulphide"] = disulphide_colour
        for i in range(len(polygons_keyslist)):
            if "X" in canvas_polygons.get(polygons_keyslist[i])[1]:
                lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        bond_colours[1] = bond_colour
        bond_colours[1] = new_colour
        for i in range(len(polygons_keyslist)):
            if canvas_polygons.get(polygons_keyslist[i])[1] == "C":
                lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        bond_colours[1] = bond_colour
        bond_colours[1] = new_colour
        for i in range(len(polygons_keyslist)):
            if canvas_polygons.get(polygons_keyslist[i])[1] == "-":
                lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        bond_colours[3] = linker_colour
        bond_colours[1] = new_colour
        for i in range(len(polygons_keyslist)):
            if canvas_polygons.get(polygons_keyslist[i])[1] == "-L-":
                lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        bond_colours[2] = hinge_colour
        bond_colours[1] = new_colour
        for i in range(len(polygons_keyslist)):
            if canvas_polygons.get(polygons_keyslist[i])[1] == "-H-":
                lower_canvas.itemconfig(polygons_keyslist[i], fill = new_colour)
        bond_colours[0] = disulphide_colour
        disulphides_list = lower_canvas.find_withtag("disulphide")
        disulphides_dict = {}
        for i in range(len(disulphides_list)):
            for j in range(len(polygons_keyslist)):
                if disulphides_list[i] == polygons_keyslist[j]:
                    disulphides_dict[j] = canvas_polygons.get(polygons_keyslist[j])
        disulphides_keyslist = list(disulphides_dict.keys())
        for i in range(len(disulphides_keyslist)):
            lower_canvas.itemconfig(disulphides_keyslist[i], fill = new_colour)
        if   selected_domain == "VHa":
            colour_checker.config(bg=a_heavy_colour)
        elif selected_domain == "VLa":
            colour_checker.config(bg=a_light_colour)
        elif selected_domain == "VHb":
            colour_checker.config(bg=b_heavy_colour)
        elif selected_domain == "VLb":
            colour_checker.config(bg=b_light_colour)
        elif selected_domain == "VHc":
            colour_checker.config(bg=c_heavy_colour)
        elif selected_domain == "VLc":
            colour_checker.config(bg=c_light_colour)
        elif selected_domain == "VHd":
            colour_checker.config(bg=d_heavy_colour)
        elif selected_domain == "VLd":
            colour_checker.config(bg=d_light_colour)
        elif selected_domain == "VHe":
            colour_checker.config(bg=e_heavy_colour)
        elif selected_domain == "VLe":
            colour_checker.config(bg=e_light_colour)
        elif selected_domain == "VHf":
            colour_checker.config(bg=f_heavy_colour)
        elif selected_domain == "VLf":
            colour_checker.config(bg=f_light_colour)
        elif selected_domain == "VHg":
            colour_checker.config(bg=g_heavy_colour)
        elif selected_domain == "VLg":
            colour_checker.config(bg=g_light_colour)
        elif selected_domain == "VHh":
            colour_checker.config(bg=h_heavy_colour)
        elif selected_domain == "VLh":
            colour_checker.config(bg=h_light_colour)
        elif selected_domain == "X":
            colour_checker.config(bg=X_colour)
        elif selected_domain == "C":
            colour_checker.config(bg=C_colour)
        elif selected_domain == "-":
            colour_checker.config(bg=bond_colour)
        elif selected_domain == "-L-":
            colour_checker.config(bg=linker_colour)
        elif selected_domain == "-H-":
            colour_checker.config(bg=hinge_colour)
        elif selected_domain == "Disulphide":
            colour_checker.config(bg=disulphide_colour)

    def colourchange_select(e):
        global specificity_colours
        global bond_colours
        global selected_colour
        global selected_domain
        global colourindex
        i=ColoursettingsLibrary.curselection()
        colourindex = i[0]
        current_colour = colours.get(coloursettings_keyslist[colourindex])
        colour_checker.config(bg=current_colour)
        selected_colour = current_colour
        selected_domain = str(coloursettings_keyslist[colourindex])





    ColoursettingsLibrary= tk.Listbox(tab2, selectbackground='#D3D3D3', height=20)
    colours = {
    "VHa": specificity_colours[0],
    "VLa":specificity_colours[1]  ,
    "VHb":specificity_colours[2] ,
    "VLb": specificity_colours[3],
    "VHc": specificity_colours[4],
    "VLc":specificity_colours[5]  ,
    "VHd":specificity_colours[6] ,
    "VLd": specificity_colours[7],
    "VHe": specificity_colours[8],
    "VLe":specificity_colours[9]  ,
    "VHf":specificity_colours[10] ,
    "VLf": specificity_colours[11],
    "VHg": specificity_colours[12],
    "VLg":specificity_colours[13]  ,
    "VHh":specificity_colours[14] ,
    "VLh": specificity_colours[15],
    "X"  : specificity_colours[18],
    "C"  : specificity_colours[19],
    "-"  : bond_colours[1],
    "-L-": bond_colours[3],
    "-H-": bond_colours[2],
    "Disulphide": bond_colours[0]
    }

    coloursettings_keyslist= list(colours.keys())
    for i in range(len(coloursettings_keyslist)):
        ColoursettingsLibrary.insert("end",coloursettings_keyslist[i])
        if i %2==0:
            ColoursettingsLibrary.itemconfig(i, bg='#D3D3D3')

    ColoursettingsLibrary.bind('<<ListboxSelect>>', colourchange_select)
    ColoursettingsLibrary.place(relx=0.01, rely = 0.1,relheight = 0.8, relwidth = 0.3)

    changerframe = tk.Frame(tab2, bg = "#D3D3D3")
    changerframe.place(relx=0.4, rely = 0.1,relheight = 0.8, relwidth = 0.49)

    colour_checker = tk.Frame(changerframe, pady = 5, bd = 5,highlightbackground="black", highlightthickness=1)
    colour_checker.place(relx=0.2, rely = 0.2,relheight = 0.2, relwidth = 0.6)

    changecolourbutton = tk.Button(changerframe, font=40, text = "Change colour", command =lambda: browse_colour())
    changecolourbutton.place(relx=0.25, rely = 0.5,relheight = 0.1, relwidth = 0.5)

    revertcolorbutton = tk.Button(changerframe, font=40, text = "Revert colour", command =lambda: revertcolor())
    revertcolorbutton.place(relx=0.25, rely = 0.7,relheight = 0.1, relwidth = 0.5)

    revertallcoloursbutton= tk.Button(changerframe, font=40, text = "Revert all colours", command =lambda: revertallcolours())
    revertallcoloursbutton.place(relx=0.25, rely = 0.9,relheight = 0.1, relwidth = 0.5)

#Coloursettings.place(relx=0.21, rely = 0.425,relheight = 0.2, relwidth = 0.19)

###Results canvas

lower_frame = tk.Frame(root, bg = '#80c1ff', bd=5)
lower_frame.place(relx=0.45, rely=0.015, relwidth=0.55,relheight=0.93)
#lower_frame.update()
frame_width = lower_frame.winfo_width()
frame_height= lower_frame.winfo_height()
print(frame_width)
print(frame_height)
lower_frame2 = tk.Frame(lower_frame, width =  700, height = 700, bg = '#80c1ff', bd=5)
lower_frame2.place(relwidth=1,relheight=1)
lower_canvas = tk.Canvas(lower_frame2,width=700,height=700+300, scrollregion=(0,0,0,700+300))
scrollbar = tk.Scrollbar(lower_frame2, orient="vertical")
scrollbar.config(command=lower_canvas.yview)
scrollbar.pack(side="right",fill="y")
lower_canvas.config(yscrollcommand=scrollbar.set)
lower_canvas.place(relheight=1,relwidth=1)

mm = MouseMover()
canvas_polygons = {}
canvas_labels   = {}
temp_label      = {}
TYPE_labels     = {}
NOTE_labels     = {}
ANTI_labels     = {}
MOD_labels      = {}
LENGTH_labels   = {}
deleted_polygons= {}
TYPE_to_redo    = {}
NOTE_to_redo    = {}
ANTI_to_redo    = {}
MOD_to_redo     = {}
LENGTH_to_redo  = {}
Polygon_to_redo = {}
Deletes_to_redo = {}
# Bind mouse events to methods (could also be in the constructor)
lower_canvas.bind("<Button-1>", mm.select)
lower_canvas.bind("<B1-Motion>", mm.drag)
lower_canvas.bind("<ButtonRelease-1>", mm.release)
lower_canvas.bind("<Button-2>", mm.change_orientation)
startcoordinates = mm.select
newcoordinates = mm.drag

domain_buttons = [InsertVHDomainButton,InsertCH1DomainButton,InsertCH2DomainButton,InsertCH3DomainButton,InsertVLDomainButton,InsertCLDomainButton,InsertCH4DomainButton,InsertXDomainButton, InsertCDomainButton, nanobody_button]
bond_buttons = [InsertBondButton,InsertLHingeButton, InsertLinkerButton,InsertDBondButton]
specificity_buttons = [a_button,b_button,c_button,d_button,e_button,f_button,g_button,h_button]
mod_buttons = [KIH_knob,KIH_hole,Positive_charge,Negative_charge,Gly_button,Mod_button]
comments_buttons = [AntiLabelButton,TypeLabelButton,NoteLabelButton,LengthLabelButton,ModLabelButton, CustomLabelButton]
delete_buttons = [InsertDelAllButton, InsertDelClickButton]
all_buttons = domain_buttons + bond_buttons + specificity_buttons + mod_buttons + comments_buttons + delete_buttons

export_frame = tk.Frame(root, bg='#FF0000')
export_frame.place(relx=0.79, rely=0.945, relwidth=0.20,relheight=0.03)
template_file_button = tk.Button(export_frame, text = "Export template file", bg = "grey", font=40, command=lambda: Get_Template_File(lower_canvas))
template_file_button.place(relx=0, rely=0, relwidth=0.5,relheight=1)
Image_file_button = tk.Button(export_frame, text = "Export PNG", bg = "grey", font=40, command=lambda: save_as_png(lower_canvas))
Image_file_button.place(relx=0.5, rely=0, relwidth=0.5,relheight=1)
img = tk.PhotoImage(file = './AbYdraw_icon.png')
root.tk.call('wm', 'iconphoto', root._w, img)

tite = root.title('abYdraw')
menubar = tk.Menu(root)
filemenu = tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Settings", command=open_settings())
filemenu.add_command(label="Open", command=lambda: browseFiles())
filemenu.add_command(label="Save", command=lambda: save_txt_file())
filemenu.add_command(label="Export template file", command=lambda: Get_Template_File(lower_canvas))
filemenu.add_command(label="Export PNG", command=lambda: save_as_png(lower_canvas))

filemenu.add_separator()

filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)

editmenu = tk.Menu(menubar, tearoff=0)
editmenu.add_command(label="Undo", command=lambda: undo())
editmenu.add_command(label="Redo", command=lambda: redo())
menubar.add_cascade(label="Edit", menu=editmenu)

helpmenu = tk.Menu(menubar, tearoff=0)
helpmenu.add_command(label="Help Index", command=donothing)
helpmenu.add_command(label="About...", command=donothing)
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)


root.mainloop()
