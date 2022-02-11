#!/usr/bin/python3
import re
import sys
import os
import tkinter as tk
from tkinter import filedialog
from tkinter import colorchooser
import tkinter.ttk as ttk
import time
import argparse
from PIL import Image



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
    ###check_for_specificities
    if "V" in str(x) and bool(re.search("\.[a-h]", str(x))) == False:
        x = re.sub("VHH","VHH.a",str(x))
        x = re.sub("VH","VH.a",str(x))
        x = re.sub("VL","VL.a",str(x))
        x = re.sub("H\.aH","HH",str(x))
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
                h_mod = 1
            else:
                h_mod = 0


            if domain == "L":
                domain = "Linker"
            note=""
            if "[" in str(chain[j]) and "]" in str(chain[j]):
                note_list= str(re.findall("\[.*?\]", str(chain[j])))
                note = str(re.sub("\[|\]|\'","", str(note_list)))

            location    = []
            locationstr =  re.findall("\((.*?)\)", str(chain[j]))
            locationstr =  str(re.sub("\[|\'|\]","", str(locationstr)))
            locationstr =  str(re.sub(":",",", str(locationstr)))
            locationstr = list(locationstr.split(","))
            if len(locationstr) == 0:
                error_message = str("ERROR: missing numbering at domain "+domain+domain+str([j])+"\nAll domains must be numbered sequentially from N-terminus to C-terminus")
                raise_error(lower_canvas, error_message)
            for i in range(len(locationstr)):
                if i == 0:
                    location_counting.append(int(locationstr[i]))
                location.append(int(locationstr[i]))





            if re.findall("\{.*?\}", str(chain[j])) != []:
                disulphide_bridges = re.findall("\{.*?\}", str(chain[j]))
                disulphide_bridges = re.sub("\{|\'|\}|\[|\]","", str(disulphide_bridges))
                disulphide_bridges = int(disulphide_bridges)
            else:
                disulphide_bridges = 0

            location = [location, disulphide_bridges,note,h_mod]
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
    if "a[" in str(VHa_keyslist+VLa_keyslist+VHb_keyslist+VLb_keyslist) and "a[" not in (fragment1_keyslist+fragment2_keyslist+fragment3_keyslist+fragment4_keyslist):
        chains = [VHa_keyslist,VLa_keyslist,VHb_keyslist,VLb_keyslist,fragment1_keyslist,fragment2_keyslist,fragment3_keyslist,fragment4_keyslist]
    else:
        chains = [VHa_keyslist,VLa_keyslist,VHb_keyslist,VLb_keyslist,fragment1_keyslist,fragment2_keyslist,fragment3_keyslist,fragment4_keyslist]

    dicts = [VHa,VLa,VHb,VLb,fragment1,fragment2,fragment3,fragment4]
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
    if ("a" not in (VHa_checked) and "a" in str(VHb_checked)):
        if "C[" in str(VHb_checked) and "C[" not in str(VHa_checked):
            pass
        else:
            VHb_checked, VHa_checked = VHa_checked,VHb_checked
            VLb_checked, VLa_checked = VLa_checked,VLb_checked
    interacting_fragments = False
    unused = []
    for i in range(len(dicts)):
        found = False
        for j in range(len(used)):
            if dicts[i] == used[j]:
                found = True
        if found == False:
            unused.append(dicts[i])
##Check claw
    Claw = False
    faux_claw = False
    Claw_number = 0
    if chain_count > 1:
        for i in range(len(chains)):
            if chains[i] != []:

                if  "X"  in chains[i][-1] or "C[" in chains[i][-1] and ("LEUCINE") not in str(dicts[i].get(chains[i][-1])[2]) and ("leucine") not in str(dicts[i].get(chains[i][-1])[2]):
                    Current_Claw_number = dicts[i].get(chains[i][-1])[0][0]
                    if Claw_number == 0:
                        Claw_number = Current_Claw_number

                    elif Claw_number != 0:
                        if Current_Claw_number == Claw_number:
                            Claw = True
                            faux_claw = True
                        if Current_Claw_number != Claw_number:
                            Claw = False
                elif  "X" not in chains[i][-1] or "C[" not in chains[i][-1]:
                    Claw = False
    #if faux_claw == True:
    #    interacting_fragments = True
    if Claw == True:
        VHa_checked = VHa
        VLa_checked = VLa
        VHb_checked = VHb
        VLb_checked = VLb
        fragment1_checked = fragment1
        fragment2_checked = fragment2
        fragment3_checked = fragment3
        fragment4_checked = fragment4
        interacting_fragments = True

    if chain_count >= 4 and Claw == False:
        if VHa_checked == {} or VHb_checked == {} or VLa_checked == {} or VLb_checked == {}:
            print(VHa_checked, VHb_checked , VLa_checked ,VLb_checked)
            error_message = "ERROR: There has been an error in pairing chains in your antibody"
            raise_error(lower_canvas, error_message)
        if fragment1 !={} and fragment2 !={} and fragment3 != {} and fragment4 !={}: #Check for second IgG
        ##Get chains that haven't been used yet
            dicts = [fragment1,fragment2,fragment3,fragment4]
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
                                            if "X" in chains[a][b] or "C[" in str(chains[a][b]):
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
        elif fragment1 !={} and fragment2 !={} and fragment3 == {} and fragment4 == {} and "X" in fragment1_keyslist[-1] and "X" not in fragment2_keyslist[-1]:
            frag1_frag2_found = False
            for i in range(len(fragment1_keyslist)):
                try:
                    current_interactor = fragment1.get(fragment1_keyslist[i])[0][1]
                    for j in range(len(fragment2_keyslist)):
                        interactor = fragment2.get(fragment2_keyslist[i])[0][0]
                        if current_interactor == interactor:
                            interacting_fragments = True
                            fragment1_checked = fragment1
                            fragment2_checked = fragment2
                            fragment3_checked = {}
                            fragment4_checked = {}
                except IndexError:
                    continue
        #elif fragment1 !={} and fragment2 !={} and fragment3 == {} and fragment4 == {} and "X" in fragment1_keyslist[-1] and "X" in fragment2_keyslist[-1]:

    all_to_check_keys = list(VHa_checked.keys())+list(VLa_checked.keys())+list(VHb_checked.keys())+list(VLb_checked.keys())+list(fragment1.keys())+list(fragment2.keys())+list(fragment3.keys())+list(fragment4.keys())
    for i in range(len(all_to_check_keys)):
        possible_domains = ["VH","VL","CH1","CH2","CH3","CH4","CL","X","H","Linker","L","C", "VHH"]
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
            domain = re.sub("\.|nano|nno|[a-h]|\@|>|\+|\-|\_|\!|\*","", domain_to_print)
            if domain not in possible_domains:
                error_message = str("ERROR: Unrecognised domain type "+ str(domain_to_print)+"\nAll domains in expression much be of type VH,VL,CH1,CH2,CH3,CH4,CL,X,H or L")
                raise_error(lower_canvas, error_message)
    if interacting_fragments == False:
        fragment1_checked = fragment1
        fragment2_checked = fragment2
        fragment3_checked = fragment3
        fragment4_checked = fragment4
    print(VHa_checked)
    print(VLa_checked)
    print(VHb_checked)
    print(VLb_checked)
    print(fragment1_checked)
    print(fragment2_checked)
    print(fragment3_checked)
    print(fragment4_checked)
    if  ((fragment1 !={} and fragment2 !={} and fragment3 != {} and fragment4 !={}) or (fragment1 !={} and fragment2 !={} and fragment3 == {} and fragment4 == {}) or Claw == True or faux_claw == True) and interacting_fragments == True:
        IgG2 = True
        if "C[" in str(fragment3_checked) and "C[" not in str(fragment1_checked):
            fragment1_checked,fragment3_checked = fragment3_checked,fragment1_checked
            fragment2_checked,fragment4_checked = fragment4_checked,fragment2_checked

    else:
        IgG2 = False
    return(VHa_checked,VLa_checked,VHb_checked,VLb_checked,chain_count,fragment1_checked,fragment2_checked,fragment3_checked,fragment4_checked,IgG2, Claw)

######################################
def Check_interactions(chains_list,canvas):
#########Set Variables################

    VHa_chain       = chains_list[0]
    VLa_chain       = chains_list[1]
    VHb_chain       = chains_list[2]
    VLb_chain       = chains_list[3]
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
    Claw            = chains_list[10]
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
    print("MULTIMERS", multimers)

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
                    ADCs.append(fragments[x][i])
                elif "C[" in fragments[x][i]:
                    CCs.append(fragments[x][i])

            fragments[x] = current_fragment

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



    def innie_or_outie(chain,VHa_chain,VHb_chain,VLa_chain,VLb_chain,Build_in,Build_out,fragment1,fragment2,fragment3,fragment4, righthanded, IgG2, Claw):
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
                                current_interactor = current_fragment.get(fragments_keyslist[j][f])[0][0]
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
                                                    elif interactor_in_or_out == "outie" and n+1!=len(keyslist):
                                                        innie_or_outie_list.append("outie")
                                                    elif interactor_in_or_out == "outie" and n+1==len(keyslist):
                                                        innie_or_outie_list.append("innie")

                                                else:
                                                    innie_or_outie_list.append("innie")
                                            else:
                                                innie_or_outie_list.append("outie")
                                        except IndexError:
                                            if "X" in keyslist[0] and VLa_chain == {} and VLb_chain == {} and "X" not in keyslist[n]:
                                                innie_or_outie_list.append("innie")
                                            else:
                                                if innie_or_outie_list[-2] == "outie":
                                                    innie_or_outie_list.append("outie")
                                                elif innie_or_outie_list[-2] == "innie":
                                                    innie_or_outie_list.append("innie")
                                                else:
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
            if "V" in keyslist[n]:
                if chain == VHa_chain or chain == VHb_chain:
                    Light_chain_check = False
                    default = "outie"
                elif chain == VLa_chain or chain == VLb_chain:
                    Light_chain_check = True
                    default = "innie"
                elif IgG2 == True and (chain == fragment1 or chain == fragment3) and Claw == False:
                    Light_chain_check = False
                    default = "outie"
                elif IgG2 == True and Claw == True:
                    Light_chain_check = True
                    Build_in = False
                    Build_out= True
                    default = "innie"
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
                    elif Claw == True:
                        in_or_out_light(n,chain,VHa_chain,Build_in,Build_out,Light_chain_check)
                    elif IgG2 == True and (chain == fragment1 or chain == fragment3):
                        innie_or_outie_list.append("outie")
                    #elif IgG2 == True and (chain == fragment2 or chain == fragment3):
                    #    if innie_or_outie_list[-2] == "outie":
                    #        innie_or_outie_list.append("innie")
                    #    elif innie_or_outie_list[-2] == "innie":
                    #        innie_or_outie_list.append("outie")
                    elif IgG2 == False and (chain == fragment1 or chain == fragment2 or chain == fragment3 or chain == fragment4) and "Linker[" in str(keyslist[n-1]):
                        if innie_or_outie_list[-2] == "outie":
                            innie_or_outie_list.append("innie")
                        elif innie_or_outie_list[-2] == "innie":
                            innie_or_outie_list.append("outie")
                    else:
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


    def renderchains(dictionary,startx,starty,canvas):
        global CLI
        global lower_canvas
        width = canvas.winfo_width()
        height = canvas.winfo_height()
        if width == 1 and height == 1 and CLI == True:
            width = 1000
            height = 900
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
        h_mods = []
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
                        if "LEUCINE" in str(dictionary) or "X" in keyslist[-1] or "X" in keyslist[0]:
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
                    try:
                        #if (dictionary == VHb_chain and "X" in keyslistb[0]) or (dictionary == VHa_chain and "X" in keyslista[0]):
                        #    Build_in = True
                        #    Build_out = False
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
                if "X" in keyslistb[0] or  "X" in keyslista[0]:
                    Build_in = True
                    Build_out = False
                else:
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
        if Claw == True:
            slant = False
            if dictionary == fragment1 or dictionary == fragment3:
                righthanded = False
            elif dictionary == fragment2 or dictionary == fragment4:
                righthanded = True


        if dictionary == VHa_chain or dictionary == VHb_chain:
            if (len(VHa_chain) == len(VHb_chain)):
                equal_chain_lengths = True
            else:
                equal_chain_lengths = False
        else:
            equal_chain_lengths = True
        innie_or_outie_list = innie_or_outie(dictionary, VHa_chain_master,VHb_chain_master,VLa_chain_master,VLb_chain_master,Build_in,Build_out, fragment1,fragment2,fragment3,fragment4,righthanded, IgG2, Claw)

        for i in range(len(dictionary)):
            keyslist = list(dictionary.keys())
            keyslist[i] = keyslist[i]
            V  = False
            X  = False
            direction = str(innie_or_outie_list[i])
            interaction = ""
            Extra_bond=False
            previous_H =False
            mod= ""
            mod_label=""
            h_mod = False
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
            if ("H[" in keyslist[i]) and dictionary.get(keyslist[i])[3] == 1:
                h_mod = True
                h_mods.append(True)
            else:
                h_mods.append(False)
            if dictionary.get(keyslist[i])[2] != "":
                note_label = str(dictionary.get(keyslist[i])[2])
                print(note_label)
                if "TYPE:" not in note_label and "NOTE:" not in note_label and "ANTI:" not in note_label and "LENGTH:" not in note_label and "MOD:" not in note_label:
                    error_message = str("ERROR: Unrecognised comment type given in ["+note_label+"]"+"\nAll comments must start with classifiers TYPE:, MOD:, NOTE:, ANTI: or LENGTH:")
                    raise_error(lower_canvas, error_message)
                if h_mod == True:
                    Domain_name = "H*"
                Notes.append(Domain_name+" "+note_label)
                if len(Notes_positions) == 0:
                    Notes_positions.append([(width/3),(height-100)])
                elif len(Notes_positions) > 0:
                    XY = (Notes_positions[-1][1])+20
                    Notes_positions.append([(width/3),XY])


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
            if CLI == False:
                print(keyslist[i])



            if i == 0:
                if "X" in keyslist[i] and dictionary.get(keyslist[i])[0] in multimers_keyslist:
                    current_number = dictionary.get(keyslist[i])[0]
                    if dictionary != VHa_chain and dictionary != VHb_chain and dictionary != VLa_chain and dictionary != VLb_chain and dictionary != fragment1 and dictionary != fragment2:
                        additions = multimers.get(dictionary.get(keyslist[i])[0])[0]
                    else:
                        additions = multimers.get(dictionary.get(keyslist[i])[0])[0]

                    if finished_multimers == []:
                        getcoordinates = domainmaker(All_positions_and_chains,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
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
                        getcoordinates = domainmaker(All_positions_and_chains,startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    else:
                        getcoordinates = domainmaker(All_positions_and_chains,startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                else:
                    if CLI == False:
                        print("checkpoint1")
                    getcoordinates = domainmaker(All_positions_and_chains,startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

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
                    if CLI == False:
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

                    getcoordinates = domainmaker(All_positions_and_chains,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                elif "H[" in keyslist[i]:
                    before_H = False
                    H_count += 1
                    Build_out = True
                    Build_in = False
                    H_coordinatey = bottom_bond

                    if i+1 ==len(dictionary):
                        if CLI == False:
                            print("checkpoint2")
                        mod = "H"
                        if dictionary == VHa_chain and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-1])[1]+3) and chain_count == 2:
                            Extra_bond=True
                            Build_up=True
                            Build_down=False
                            if slant == True and righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif slant == True and righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif slant==False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)



                        else:# or dictionary == VHb_chain:
                            if CLI == False:
                                print("checkpoint3")


                            if righthanded == True and slant==True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == False and slant==True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif slant==False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)


                elif "X" in keyslist[i] and dictionary.get(keyslist[i])[0] in multimers_keyslist:
                    current_number = dictionary.get(keyslist[i])[0]
                    additions = multimers.get(dictionary.get(keyslist[i])[0])[0]

                    if finished_multimers != []:
                        if righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(finished_multimers[-1][0]-additions[0]),(finished_multimers[-1][1]+additions[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(finished_multimers[-1][0]+additions[0]),(finished_multimers[-1][1]+additions[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    elif finished_multimers == []:
                        if righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-80),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+80),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                elif "H[" in keyslist[i-1]  and ("X" in keyslist[i] or "C[" in keyslist[i]):
                    if CLI == False:
                        print("checkpoint4")
                    if righthanded == False:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    elif righthanded == True:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)



                elif ("X[" in keyslist[i] or "C[" in keyslist[i]) and "H[" not in keyslist[i-1] and len(dictionary.get(keyslist[i])) == 1:
                    if CLI == False:
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
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True and slant == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif slant == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                    elif chain_count !=2 and mod !="Leucine" and "X" not in keyslist[i-1]:
                        if innie_or_outie_list[i-2] == "innie" and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif innie_or_outie_list[i-2] == "outie" and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif innie_or_outie_list[i-2] == "innie" and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif innie_or_outie_list[i-2] == "outie" and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                        elif righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        Build_in = True
                        Build_out = False
                    elif chain_count !=2 and mod !="Leucine" and "X" in keyslist[i-1]:
                        Build_up = True
                        if righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                        Build_in = True
                        Build_out = False
                    elif mod =="Leucine":
                        if righthanded == False and slant== True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == True and slant==True :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif slant == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)


                elif "H[" in keyslist[i-1]  and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[1] != (dictionary.get(previous_domain)[0]):
                    if CLI == False:
                        print("checkpoint6")
                    previous_H = True
                    slant=False
                    if chain_count <=2:
                        if dictionary.get(keyslist[2]) != [''] and "Linker[" in keyslist[2]:
                            try:
                                if dictionary.get(keyslist[0])[0] != (dictionary.get(keyslist[2])[1]) :
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            except IndexError:
                                if righthanded == True :
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                elif righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                        else:
                            if righthanded == True :
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                    elif H_count ==1 :
                        if righthanded == True :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    elif H_count==2 :
                        if righthanded == True :
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                    Build_in = False
                    Build_out = True

                elif "Linker[" not in keyslist[i-1] and "X" and len(dictionary.get(keyslist[i-1])) ==1:
                    if CLI == False:
                        print("checkpoint7")
                    if "Linker[" in keyslist[i]:
                        pass
                        if CLI == False:
                            print("checkpoint7.5")
                    elif slant == True and righthanded == True:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    else:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    Build_in  = False
                    Build_out = True

                elif chain_count == 2 and "H[" not in keyslist[i] and "Linker[" not in keyslist[i-1] and "Linker[" not in keyslist[i] and "X" not in keyslist[i] and "X" not in keyslist[i-1] and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-1])[1]+3):
                    if CLI == False:
                        print("checkpoint8")
                    if chain_count == 2:
                        if dictionary == VHa_chain:
                            Build_up=True
                            if slant == True and righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif slant == True and righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            else:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            Build_up=True
                            Build_down=False
                        elif dictionary == VHb_chain:
                                if slant == True and righthanded == True:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-50,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                elif slant == True and righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+50,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                else:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                    #Build_in  = False
                    #Build_out = True


                elif "H[" not in keyslist[i] and "Linker[" not in keyslist[i-1] and "Linker[" not in keyslist[i] and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):

                    if CLI == False:
                        print("checkpoint9")
                    if slant == True and righthanded == True:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                    else:
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)


                    #Build_in  = False
                    #Build_out = True

##Linker
                elif "Linker[" in keyslist[i-1]:
                    if CLI == False:
                        print("checkpoint10")
                    previous_chain = chain[i-2]
                    previous_domain = keyslist[i-2]
                    previous_number = previous_number+1


##SdFV
                    if len(dictionary.get(keyslist[i])) == 1:
                        if CLI == False:
                            print("checkpoint11")
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        else:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        Build_in  = False
                        Build_out = True
                    elif len(dictionary.get(keyslist[i-2])) ==1 and "X" not in keyslist[i-2]:
                        if CLI == False:
                            print("checkpoint11.5")
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        else:
                            getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        Build_in  = False
                        Build_out = True

##Self-Interacting chains

                    elif dictionary.get(keyslist[i])[0] == previous_number and str(dictionary.get(keyslist[i])[1]) in Location_Text:
                        if CLI == False:
                            print("checkpoint12")
                        to_join_number      = str(dictionary.get(keyslist[i])[1])
                        to_join_coordinates = find_the_fragment(to_join_number,All_positions_and_chains)
                        to_joinx            = to_join_coordinates[0][0]
                        to_joiny            = to_join_coordinates[0][1]
                        to_join_righthanded = to_join_coordinates[1]
                        to_join_direction   = to_join_coordinates[2]
                        all_list = list(All_positions_and_chains.keys())
                        get = (All_positions_and_chains.get(all_list[-1]))
                        yprev = get[0][1]
                        if dictionary.get(keyslist[i-2])[0] != dictionary.get(keyslist[i])[1] and to_joiny < yprev:
                            if CLI == False:
                                print("checkpoint13")
                            Build_up=True
                            Build_down=False
                        if chain_count == 1 and  dictionary.get(keyslist[i])[1] != (dictionary.get(keyslist[i-2])[0]) and len(dictionary) >=8 :#and "X" not in keyslist[i] :#and to_join_direction == ("innie" or "outie"):
                            print("checkpoint12.6")
                            righthanded = True
                            Build_in = False
                            Build_out = True
                        elif chain_count == 1  and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-2])[0]):
                            print("checkpoint12.7")
                            Build_in = False
                            Build_out = True
                        if to_join_righthanded == False and to_join_direction == 'outie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif to_join_righthanded == False and to_join_direction == 'innie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif to_join_righthanded == False and to_join_direction == 'constant':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif to_join_righthanded == True and to_join_direction == 'outie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif to_join_righthanded == True and to_join_direction == 'innie':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif to_join_righthanded == True and to_join_direction == 'constant':
                            getcoordinates = domainmaker(All_positions_and_chains,to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)




                        if Build_in == True:
                            if CLI == False:
                                print("checkpoint13.1")
                            Build_out =True
                            Build_in = False
                        elif Build_out==True:
                            if CLI == False:
                                print("checkpoint13.2")
                            Build_in = True
                            Build_out = False
                        #if change_side == True:
                        #    Build_out =True
                        #    Build_in = False


##ADCs
                    elif ("X" in keyslist[i-2] or "C[" in keyslist[i-2]) :
                        if CLI == False:
                            print("checkpoint14")
                        if chain_count ==1:
                            if Build_in == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                elif righthanded == True:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif Build_out == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                elif righthanded == True:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif chain_count > 1:
                            if chain_count == 2 and keyslist[i-2] != keyslist[0]:
                                slant = False
                            if righthanded == False and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == False and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == True and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == True and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            if chain_count==2 and keyslist[i-2] != keyslist[0]:
                                Build_out = True
                                Build_in = False
                            else:
                                Build_in = True
                                Build_out = False
##Build up
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and "X" not in keyslist[0] and dictionary.get(keyslist[i])[0] == (dictionary.get(keyslist[0])[1]):
                        if CLI == False:
                            print("checkpoint15")
                        Build_up=True
                        Build_down=False
                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                    elif chain_count == 2 and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and  dictionary.get(keyslist[i])[1]+2 == dictionary.get(keyslist[i-2])[1]:
                        if CLI == False:
                            print("checkpoint16")
                        if dictionary == VHa_chain:
                            Build_up=True
                            Build_down=False
                            if slant==True and righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])-45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif  slant==True and righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])+45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif slant==False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                        elif dictionary == VHb_chain:
                            try:
                                if "Linker[" in keyslist[i+1]:

                                    if slant==True and righthanded == True:
                                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])+45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                        Build_up=True
                                        Build_down=False
                                    elif  slant==True and righthanded == True:
                                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0])-45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                                        Build_up=True
                                        Build_down=False
                                    elif slant==False:
                                        getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                                else:
                                    getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            except IndexError:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

                        if Build_out == True:
                            Build_in = True
                            Build_out = False
                        elif Build_in == True:
                            Build_out = True
                            Build_in = False

##Build across

                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] == (dictionary.get(previous_domain)[1]):
                        if CLI == False:
                            print("checkpoint17")
                        if Build_in == True:
                            if righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif Build_out == True:
                            if righthanded == False:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                            elif righthanded == True:
                                getcoordinates = domainmaker(All_positions_and_chains,(previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

##Build down
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                        if CLI == False:
                            print("checkpoint18")
                        if "V" in keyslist[i]:
                            in_out_counter +=1

                        def check_for_clashes(startx, starty):
                            clash = False
                            All_positions_and_chains_list = list(All_positions_and_chains.keys())
                            for i in range(len(All_positions_and_chains_list)):
                                coordinates = All_positions_and_chains.get(All_positions_and_chains_list[i])[0]
                                if len(coordinates) > 2:
                                    x1 = get_min_max_coordinates(coordinates)[0]
                                    x2 = get_min_max_coordinates(coordinates)[1]
                                    y1 = get_min_max_coordinates(coordinates)[2]
                                    y2 = get_min_max_coordinates(coordinates)[3]
                                    if startx == coordinates[0] and starty == y1:
                                        clash = True
                            if clash == True:
                                return(True,startx,starty+95 )
                            elif clash == False:
                                return(False,startx,starty )

                        clash = check_for_clashes((previous_chain[6]),(previous_chain[7]+20))
                        if clash[0] == True:
                            Build_up = True
                            Build_down=False
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker(All_positions_and_chains,clash[1],clash[2],righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker(All_positions_and_chains,clash[1],clash[2],righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)
                        else:
                            getcoordinates = domainmaker(All_positions_and_chains,clash[1],clash[2],righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up)

##H disulphides
            if "H[" in keyslist[i] or "H*[" in keyslist[i] or "Linker[" in keyslist[i] and len(dictionary.get(keyslist[i])) >1 and (dictionary != VHa_1_test or dictionary != VHb_1_test or dictionary != VLa_1_test or dictionary != VLb_1_test):
                print("DISULPHIDE CHECKPOINT")
                if Build_up == False:
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
                        #slant = True
                        if slant == True and righthanded == False:
                            topx = bottomx-20
                        elif slant == True and righthanded == True:
                            topx = bottomx+20
                        elif slant==False:
                            topx = bottomx
                        topy = bottomy+40
                    location   = dictionary.get(keyslist[i])[0]
                    if "H" in keyslist[i]:
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
                    elif "Linker[" in keyslist[i]:
                        bottomy = bottomy-15
                        try:
                            if dictionary.get(keyslist[i-1])[0] == dictionary.get(keyslist[i+1])[0][1]:
                                bottomy += 20
                                topy = bottomy-80
                                if Build_in == True:
                                    if righthanded == False:
                                        topx = bottomx+20
                                    elif righthanded == True:
                                        topx = bottomx-20
                                elif Build_in == False:
                                    if righthanded == False:
                                        topx = bottomx-60
                                    elif righthanded == True:
                                        topx = bottomx+60
                        except IndexError:
                            try:
                                if dictionary.get(keyslist[i-1])[1] == dictionary.get(keyslist[i+1])[0][0]:
                                    bottomy += 20
                                    topy = bottomy-80
                                    if Build_in == True:
                                        if righthanded == False:
                                            topx = bottomx+20
                                        elif righthanded == True:
                                            topx = bottomx-20
                                    elif Build_in == False:
                                        if righthanded == False:
                                            topx = bottomx-60
                                        elif righthanded == True:
                                            topx = bottomx+60
                            except IndexError:
                                pass
                        location   = dictionary.get(keyslist[i])[0]
                        interactor = dictionary.get(keyslist[i])[1]
                        if bottomy > topy:
                            bottomy,topy = topy,bottomy
                            bottomx,topx = topx,bottomx
                        H_bond_coordinates = disulphide_maker(disulphide_bridge_count,bottomx,bottomy,topx,topy,righthanded)
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
                if CLI == False:
                    print("checkpoint19")
                top_bond = getcoordinates[2]
                if "H[" in keyslist[i-1]:
                    if h_mods[-2] == True:
                        h_mod= True
                    hinges.append([[bottom_bond + top_bond],h_mod])
                elif "Linker[" in keyslist[i-1]:
                    linkers.append(bottom_bond + top_bond)
                    global L_Labels
                    if L_Labels == True:
                        L_Domain_name = "L"
                        if top_bond[0] >= bottom_bond[0]:
                            L_Label_x = ((top_bond[0]+bottom_bond[0])/2)
                        elif bottom_bond[0] >= top_bond[0]:
                            L_Label_x = ((bottom_bond[0]+top_bond[0])/2)
                        if top_bond[1] >= bottom_bond[1]:
                            L_Label_y = ((top_bond[1]+bottom_bond[1])/2)
                        elif bottom_bond[1] >= top_bond[1]:
                            L_Label_y = ((bottom_bond[1]+top_bond[1])/2)
                        L_Label_Locations = [L_Label_x,L_Label_y]
                        if dictionary.get(keyslist[i])[0] == (dictionary.get(keyslist[i-2])[1]):
                            L_Label_Locations[1] -= 30
                        elif dictionary.get(keyslist[i])[0] != (dictionary.get(keyslist[i-2])[1]):
                            if righthanded == False and slant == True:
                                L_Label_Locations[0] -= 40
                            elif righthanded == False and slant == False:
                                L_Label_Locations[0] -= 20
                            elif righthanded == True and slant == True:
                                L_Label_Locations[0] += 40
                            elif righthanded == True and slant == False:
                                L_Label_Locations[0] += 20
                        L_text = "L"
                        text_coordinates.append([L_Label_Locations])
                        Location_Text.append(str(L_text)+mod_label)
                        Domain_Text.append(str(L_Domain_name)+mod_label)

                else:

                    bonds.append(bottom_bond + top_bond)
                #if mod=="Leucine":
                #    bonds.append(getcoordinates[0])
            elif i > 0 and Build_up==True:
                if CLI == False:
                    print("checkpoint20")
                changed_righthand = False
                if chain_count == 1 and righthanded == True:
                    righthanded = False
                    changed_righthand = True
                top_bond = getcoordinates[2]
                arc_topx  = top_bond[0]
                arcbottomx= bottom_bond[0]
                if arc_topx != arcbottomx and slant == False:
                    top_bond = getcoordinates[2]
                    arc_topx  = top_bond[0]
                    arc_topy  = top_bond[1]
                    arcbottomx= bottom_bond[0]
                    arcbottomy= bottom_bond[1]
                    if "Linker[" not in keyslist[i-1]:
                        if righthanded == True:
                            bonds.append([bottom_bond[0], bottom_bond[1], top_bond[0]-20,top_bond[1]])
                        elif righthanded == False:
                            bonds.append([bottom_bond[0], bottom_bond[1], top_bond[0]+20,top_bond[1]])
                    elif "Linker[" in keyslist[i-1]:
                        if righthanded == True:
                            linkers.append([bottom_bond[0], bottom_bond[1], top_bond[0]+3,top_bond[1]-18])
                        elif righthanded == False:
                            linkers.append([bottom_bond[0], bottom_bond[1], top_bond[0]-3,top_bond[1]-18])

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
                            hinges.append([[extra_bond], h_mod])
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



                if changed_righthand == True:
                    righthanded = True





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
                if CLI == False:
                    print("checkpoint27")
                Label_Locations = getcoordinates[3]
                text_coordinates.append([Label_Locations])
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [getcoordinates[0], righthanded,direction]
                Domain_Text.append(str(Domain_name)+mod_label)


            elif ("X[" in keyslist[i] or "C[" in keyslist[i]):
                global Show_Leucine_Zippers
                if  "X[" in keyslist[i] and (mod != "Leucine" or Show_Leucine_Zippers == False):
                    ADCs.append(getcoordinates[0])
                elif "C[" in keyslist[i] and (mod != "Leucine" or Show_Leucine_Zippers == False):
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

            if dictionary != VHa_chain and dictionary != VHb_chain and dictionary != VLa_chain and dictionary != VLb_chain and dictionary != fragment1 and dictionary != fragment2 and dictionary != fragment3 and dictionary != fragment4 and "Linker[" not in keyslist[i]:

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

                    if "Linker[" not in keyslist[i]:
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
                        hinges.append([[bottomx,bottomy,topx,topy],h_mod])
                    elif slant == False:
                        if tangle_found == True and dictionary == VHa_chain:
                            pass
                        else:
                            hinges.append([[bottomx,bottomy,topx,topy+20],h_mod])

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
                text = dictionary.get(keyslist[i])[0]
                All_positions_and_chains[text] = [Label_Locations, righthanded, innie_or_outie, direction]
                global H_Labels
                if H_Labels == True:
                    text_coordinates.append(Label_Locations)
                    Location_Text.append(str(text)+mod_label)
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
    global CLI
    width = canvas.winfo_width()
    height = canvas.winfo_height()
    if width == 1 and height == 1 and CLI == True:
        width = 1000
        height = 900
    if chain_count == 1:
        VHa_1_test = VHa_chain.copy()
        VLa_1_test = VLa_chain.copy()
        VHb_1_test = VHb_chain.copy()
        VLb_1_test = VLb_chain.copy()
        VHa_startx, VHa_starty = (width/2),(height/2)-200
        VHb_startx, VHb_starty = 0,0
        VLa_startx, VLa_starty = 0,0
        VLb_startx, VLb_starty = 0,0
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty, canvas)
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty, canvas)

    elif chain_count >= 3 or (chain_count == 2  and tangle_found == False and ("H[" in str(VHa_chain) and "H[" in str(VHb_chain)) or ("X[" in str(VHa_chain) and "X[" in str(VHb_chain)) or ("C[" in str(VHa_chain) and "C[" in str(VHb_chain))):
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
                            if Xa == Xb or Claw == True:
                                VHa_H_coordinatesx = (width/2)
                                VHa_H_coordinatesy = (height/2)-50
                                VHb_H_coordinatesx = (width/2)
                                VHb_H_coordinatesy = (height/2)-50
                            elif Xa != Xb:
                                VHa_H_coordinatesx = (width/2)-30
                                VHa_H_coordinatesy = (height/2)-50
                                VHb_H_coordinatesx = (width/2)+30
                                VHb_H_coordinatesy = (height/2)-50

        elif chain_count == 4 and "H[" not in str(VHa_1_test) and "H[" not in str(VHb_1_test) and "X[" not in str(VHa_1_test) and "X[" not in str(VHb_1_test) and "C[" not in str(VHa_1_test) and "C[" not in str(VHb_1_test):
            VHa_H_coordinatesx = (width/2)-32
            VHa_H_coordinatesy = (height/2)
            VHb_H_coordinatesx = (width/2)+32
            VHb_H_coordinatesy = (height/2)
        teststartx = 0
        teststarty = 0
        testHpositionVHa = renderchains(VHa_1_test,teststartx,teststarty, canvas)[25]
        test_H_positionx = testHpositionVHa[0]
        test_H_positiony = testHpositionVHa[1]
        differencetest_desiredx = VHa_H_coordinatesx - test_H_positionx
        differencetest_desiredy = VHa_H_coordinatesy - test_H_positiony
        VHa_startx = teststartx + differencetest_desiredx
        VHa_starty = teststarty + differencetest_desiredy

##VHb_chain

        teststartx = width
        teststarty = 0
        testHpositionVHb = renderchains(VHb_1_test,teststartx,teststarty, canvas)[25]
        test_H_positionx = testHpositionVHb[0]
        test_H_positiony = testHpositionVHb[1]
        differencetest_desiredx1 = test_H_positionx - VHb_H_coordinatesx
        differencetest_desiredx2 =  VHb_H_coordinatesx - test_H_positionx
        if differencetest_desiredx1 > differencetest_desiredx2:
            differencetest_desiredx = differencetest_desiredx1
        elif differencetest_desiredx1 < differencetest_desiredx2:
            differencetest_desiredx = differencetest_desiredx2
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
            #VHb_startx -= 515
            VHb_startx = (width/5)*3.3
        elif VHb_startx > width and IgG2 == False:
            VHb_startx -= 600
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty, canvas)
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty, canvas)


###Get start positions of light chains and render
    elif chain_count == 2 and Claw == False:

        keyslista = list(VHa_chain.keys())
        keyslistb = list(VHb_chain.keys())
        keyslist = list(VHa_chain.keys())
        VHb_startx, VHb_starty = (width/2)+100,(height/2)-200
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty, canvas)
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

        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty, canvas)



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
                        testHpositionVLa = renderchains(VLa_1_test,teststartx,teststarty, canvas)[25]
                        test_H_positionx = testHpositionVLa[0]
                        test_H_positiony = testHpositionVLa[1]
                        differencetest_desiredx = VLa_refx - test_H_positionx
                        differencetest_desiredy = VLa_refy - test_H_positiony
                        VLa_startx = teststartx + differencetest_desiredx
                        VLa_starty = teststarty + differencetest_desiredy
                        break
        if interactor_count == 0 and Claw == True:
            VLa_startx, VLa_starty = VHa_startx,VHa_starty+100
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
                        testHpositionVLb = renderchains(VLb_1_test,teststartx,teststarty, canvas)[25]
                        test_H_positionx = testHpositionVLb[0]
                        test_H_positiony = testHpositionVLb[1]
                        differencetest_desiredx = VLb_refx - test_H_positionx
                        differencetest_desiredy = VLb_refy - test_H_positiony
                        VLb_startx = teststartx + differencetest_desiredx
                        VLb_starty = teststarty + differencetest_desiredy
                        break
        if interactor_count == 0 and Claw == True:
            VLb_startx,VLb_starty = VHb_startx,VHb_starty+100

    else:
        VLb_startx,VLb_starty = 0,0

    VLa_stats = renderchains(VLa_chain,VLa_startx, VLa_starty, canvas)
    VLb_stats = renderchains(VLb_chain,VLb_startx, VLb_starty, canvas)


##Two-chained Abs


    frag1_startx,frag1_starty,frag2_startx,frag2_starty,frag3_startx,frag3_starty,frag4_startx,frag4_starty = 0,0,0,0,0,0,0,0
    if IgG2 == False:
        def place_fragments(fragment, backup_startx,backup_starty):
            All_positions_and_chains_keys = list(All_positions_and_chains.keys())
            if fragment != {}:
                fragment_list = list(fragment.keys())
                interaction_found = False
                connection_found = False

                for i in range(len(fragment_list)):
                    try:
                        fragment_inter = fragment.get(fragment_list[i])[0][1]
                        frag_start = find_the_fragment(fragment_inter,All_positions_and_chains)
                        if frag_start is not None and interaction_found == False:
                            interaction_found = True
                            righthanded = frag_start[1]
                            break
                    except IndexError:
                        inter = fragment.get(fragment_list[i])[0][0]
                        for i in range(len(All_positions_and_chains_keys)):
                            if All_positions_and_chains_keys[i] == inter:
                                connection_found = True

                if interaction_found == True:
                    if righthanded == True:
                        frag_startx=frag_start[0][0]+60
                        frag_starty=frag_start[0][1]
                    elif righthanded==False:
                        frag_startx=frag_start[0][0]-60
                        frag_starty=frag_start[0][1]

                elif interaction_found == False and connection_found == True:
                    frag_startx, frag_starty = backup_startx,backup_starty

                elif interaction_found == False and connection_found == False:
                    error_message = "ERROR:Your expression contains multiple antibodies that are not connected"
                    raise_error(lower_canvas,error_message)

                return(frag_startx, frag_starty)
            else:
                return(0,0)

        frag1_xy = place_fragments(fragment1, VHa_startx-100,VHa_starty+200)
        frag2_xy = place_fragments(fragment2, VLa_startx-100,VHa_starty+350)
        frag3_xy = place_fragments(fragment3, VHb_startx-100,VHb_starty+200)
        frag4_xy = place_fragments(fragment4, VLb_startx-100,VLb_starty+350)
        frag1_stat= renderchains(fragment1,frag1_xy[0],frag1_xy[1], canvas)
        frag2_stat= renderchains(fragment2,frag2_xy[0],frag2_xy[1], canvas)
        frag3_stat= renderchains(fragment3,frag3_xy[0],frag3_xy[1], canvas)
        frag4_stat= renderchains(fragment4,frag4_xy[0],frag4_xy[1], canvas)
    elif IgG2 == True:
        if Claw == True:
            frag1_stat= renderchains(fragment1,VHa_startx-0,VHa_starty+200, canvas)
            frag2_stat= renderchains(fragment2,VHb_startx-0,VHb_starty+200, canvas)
            frag3_stat= renderchains(fragment3,VHa_startx-0,VHa_starty+300, canvas)
            frag4_stat= renderchains(fragment4,VHb_startx-0,VHb_starty+300, canvas)

        elif ((("X" in str(list(fragment1.keys())) and "X" in str(list(fragment3.keys()))) or ("C[" in str(list(fragment1.keys())) and "C[" in str(list(fragment3.keys())))) and ("CH2" not in str(list(fragment1.keys())) and "CH2" not in str(list(fragment3.keys())))) or fragment3 == {} and fragment4 == {} :
            frag1_stat= renderchains(fragment1,VHa_startx-0,VHa_starty+200, canvas)
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
            try:
                righthanded = frag2_start[1]
                if righthanded == True:
                    frag2_startx=frag2_start[0][0]+60
                    frag2_starty=frag2_start[0][1]
                elif righthanded==False:
                    frag2_startx=frag2_start[0][0]-60
                    frag2_starty=frag2_start[0][1]
                frag2_stat = renderchains(fragment2,frag2_startx,frag2_starty, canvas)
            except TypeError:
                frag2_stat = renderchains(fragment2,frag2_startx,frag2_starty, canvas)
            if fragment3 != {} and fragment4 != {}:
                frag3_stat= renderchains(fragment3,VHb_startx+0,VHb_starty+200, canvas)
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
                frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty, canvas)
            else:
                frag3_stat = renderchains(fragment3,0,0, canvas)
                frag4_stat= renderchains(fragment4,0,0, canvas)

        elif fragment1 != {} and fragment2 !={} and fragment3 !={} and fragment4 !={}:

            VHa_H_coordinatesx = ((width/6)*4)-10
            VHa_H_coordinatesy = (height/2)+100
            VHb_H_coordinatesx = VHa_H_coordinatesx+100
            VHb_H_coordinatesy = (height/2)+100

            frag1_stat= renderchains(fragment1,VHa_startx+325,VHa_starty+200, canvas)
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
            conj_bond = []
            for i in range(len(coordinates_to_change)):
                for j in range(len(coordinates_to_change[i])):
                    for k in range(len(coordinates_to_change[i][j])):
                        if coordinates_to_change[i][j] != conj_bond:
                            if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
                                try:
                                    if len(coordinates_to_change[i][j]) == 4 and coordinates_to_change[i][j][1] > conj_y2 and conj_x1 <= coordinates_to_change[i][j][2] <= conj_x2 and conj_y1 <= coordinates_to_change[i][j][3] <= conj_y2 and conj_fixed == False:
                                        if frag1_differencetest_desiredy <= 0 and frag1_differencetest_desiredx >= 0:
                                            coordinates_to_change[i][j][0]+= frag1_differencetest_desiredx
                                            coordinates_to_change[i][j][1]+= frag1_differencetest_desiredy
                                        elif frag1_differencetest_desiredy >= 0 and frag1_differencetest_desiredx <= 0:
                                            coordinates_to_change[i][j][0]-= frag1_differencetest_desiredx
                                            coordinates_to_change[i][j][1]-= frag1_differencetest_desiredy
                                        elif frag1_differencetest_desiredy <= 0 and frag1_differencetest_desiredx <= 0:
                                            coordinates_to_change[i][j][0]-= frag1_differencetest_desiredx
                                            coordinates_to_change[i][j][1]+= frag1_differencetest_desiredy
                                        elif  frag1_differencetest_desiredy >= 0 and frag1_differencetest_desiredx >= 0:
                                            coordinates_to_change[i][j][0]+= frag1_differencetest_desiredx
                                            coordinates_to_change[i][j][1]+= frag1_differencetest_desiredy
                                        if "C[" in str(VHa_chain) and "C[" not in str(VHb_chain):
                                            coordinates_to_change[i][j][2] = VHa_stats[52][0][0]
                                            coordinates_to_change[i][j][3] = VHa_stats[52][0][1]+10
                                        elif "C[" in str(VHb_chain) and "C[" not in str(VHa_chain):
                                            coordinates_to_change[i][j][2] = VHb_stats[52][0][0]
                                            coordinates_to_change[i][j][3] = VHb_stats[52][0][1]+10
                                        conj_fixed = True
                                        conj_bond = coordinates_to_change[i][j]

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

            #print(frag1_stat)
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
            frag2_stat = renderchains(fragment2,frag2_startx,frag2_starty, canvas)

            frag3_stat= renderchains(fragment3,VHb_startx+325,VHb_starty+200, canvas)
            test_H_positionfrag3 = frag3_stat[26]


            test_H_positionx = test_H_positionfrag3[0][0]
            test_H_positiony = test_H_positionfrag3[0][1]
            frag3_differencetest_desiredx = VHb_H_coordinatesx- test_H_positionx
            frag3_differencetest_desiredy = VHb_H_coordinatesy- test_H_positiony
            coordinates_to_change = [frag3_stat[0],frag3_stat[1],frag3_stat[2],frag3_stat[3],frag3_stat[4],frag3_stat[5],frag3_stat[6],frag3_stat[7],frag3_stat[36],frag3_stat[37],frag3_stat[38],frag3_stat[39],frag3_stat[40],frag3_stat[41],frag3_stat[42],frag3_stat[43],frag3_stat[8],frag3_stat[26],frag3_stat[27],frag3_stat[24],frag3_stat[11], frag3_stat[10]]
            for i in range(len(coordinates_to_change)):
                for j in range(len(coordinates_to_change[i])):
                    for k in range(len(coordinates_to_change[i][j])):
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
            frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty, canvas)
            hingey1 = frag3_stat[26][0][1]


            for i in range(len(completed_disulphidebridges)):
                if completed_disulphidebridges[i][1] > hingey1:
                    completed_disulphidebridges[i][2] -= frag1_differencetest_desiredx
                    completed_disulphidebridges[i][3] -= frag1_differencetest_desiredy
                    completed_disulphidebridges[i][0] += frag3_differencetest_desiredx
                    completed_disulphidebridges[i][1] += frag3_differencetest_desiredy


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
                if All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j] < 0 :
                    Xout_of_range = True
                    if Xhow_much > All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j]:
                        Xhow_much = All_positions_and_chains.get(keyslist_All_positions_and_chains[i])[0][j]

    if Yout_of_range == True:
        Ynew_start = Yhow_much+5
    if Xout_of_range == True:
        Xnew_start = Xhow_much
    #print("X out of range", Xnew_start)




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
        Ynew_start = Yhow_much-10
        coordinates_to_change = [Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Heavy_Domains_e,Light_Domains_e,Heavy_Domains_f,Light_Domains_f,Heavy_Domains_g,Light_Domains_g,Heavy_Domains_h,Light_Domains_h,Bonds,Hinges,Linkers,ADCs,disulphide_bridges, Label_spot,Chem_con]
        for i in range(len(coordinates_to_change)):
            for j in range(len(coordinates_to_change[i])):
                for k in range(len(coordinates_to_change[i][j])):
                    if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
                        if  k %2 != 0:
                            coordinates_to_change[i][j][k] -= Ynew_start
                    else:
                        for l in range(len(coordinates_to_change[i][j][k])):
                            if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
                                if l %2 != 0:
                                    coordinates_to_change[i][j][k][l] -= Ynew_start
    #if Xout_of_range == True:
    #    Xnew_start = Xhow_much
    #    coordinates_to_change = [Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Heavy_Domains_e,Light_Domains_e,Heavy_Domains_f,Light_Domains_f,Heavy_Domains_g,Light_Domains_g,Heavy_Domains_h,Light_Domains_h,Bonds,Hinges,Linkers,ADCs,disulphide_bridges, Label_spot,Chem_con]
    #    for i in range(len(coordinates_to_change)):
    #        for j in range(len(coordinates_to_change[i])):
    #            for k in range(len(coordinates_to_change[i][j])):
    #                if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
    #                    if  k %2 == 0:
    #                        coordinates_to_change[i][j][k] -= Xnew_start
    #                else:
    #                    for l in range(len(coordinates_to_change[i][j][k])):
    #                        if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
    #                            if l %2 == 0:
    #                                coordinates_to_change[i][j][k][l] -= Xnew_start
    print(Hinges,Heavy_Domains_a,names_list_heavy_a,Light_Domains_a,names_list_light_a,Heavy_Domains_b,names_list_heavy_b,Light_Domains_b,names_list_light_b,Heavy_Domains_c,names_list_heavy_c,Light_Domains_c,names_list_light_c,Heavy_Domains_d,names_list_heavy_d,Light_Domains_d,names_list_light_d,Label_Text,Label_spot,Domain_Text,Notes,Notes_positions,arcs_left,arcs_right,arcs_left_slant,arcs_right_slant,ADCs,Heavy_Domains_e,names_list_heavy_e,Light_Domains_e,names_list_light_e,Heavy_Domains_f,names_list_heavy_f,Light_Domains_f,names_list_light_f,Heavy_Domains_g,names_list_heavy_g,Light_Domains_g,names_list_light_g,Heavy_Domains_h,names_list_heavy_h,Light_Domains_h,names_list_light_h)
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
        global Bond_thickness
        global Bond_Arrows
        Arrow_dimensions = Bond_Arrows

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
                domain = canvas.create_line(disulphide_bridges[i], fill='#FF4040', width = Bond_thickness,tags="disulphide", arrow=tk.BOTH, arrowshape=Arrow_dimensions)
                canvas_polygons[domain] = [disulphide_bridges[i], "-disulphide-"]

    #Bonds
        for i in range(len(Bonds)):
            domain = canvas.create_line(Bonds[i], fill=bond_colour, width = Bond_thickness,tags="bonds", arrow=tk.LAST, arrowshape=Arrow_dimensions)
            canvas_polygons[domain] = [Bonds[i], "-"]
        if arcs_left!=[]:
            for i in range(len(arcs_left)):
                if arcs_left[i][1] == False:
                    domain = canvas.create_arc(arcs_left[i][0], start=90, extent=180, style=tk.ARC, fill=bond_colour, width = Bond_thickness,tags=("bonds","arcs_left","arcs"))
                    canvas_polygons[domain] = [arcs_left[i][0], "-"]
                elif arcs_left[i][1] == True:
                    domain = canvas.create_arc(arcs_left[i][0], start=90, extent=180, style=tk.ARC, outline=linker_colour, width = Bond_thickness,tags=("bonds","arcs_left","arcs"))
                    canvas_polygons[domain] = [arcs_left[i][0], "-L-"]
        if arcs_right!=[]:
            for i in range(len(arcs_right)):
                if arcs_right[i][1] == False:
                    domain = canvas.create_arc(arcs_right[i][0], start=270, extent=180, style=tk.ARC, fill=bond_colour, width = Bond_thickness,tags=("bonds","arcs_right","arcs"))
                    canvas_polygons[domain] = [arcs_right[i][0], "-"]
                elif arcs_right[i][1] == True:
                    domain = canvas.create_arc(arcs_right[i][0], start=270, extent=180, style=tk.ARC, outline=linker_colour, width = Bond_thickness,tags=("bonds","arcs_right","arcs"))
                    canvas_polygons[domain] = [arcs_right[i][0], "-L-"]
        if arcs_left_slant != []:
            for i in range(len(arcs_left_slant)):
                if arcs_left_slant[i][1] == False:
                    domain = canvas.create_arc(arcs_left_slant[i][0], start=150, extent=120, fill = bond_colour, style=tk.ARC,width=Bond_thickness,tags=("bonds","arcs_left_slant","arcs"))
                    canvas_polygons[domain] = [arcs_left_slant[i][0], "-"]
                elif arcs_left_slant[i][1] == True:
                    domain = canvas.create_arc(arcs_left_slant[i][0], start=150, extent=120, outline = linker_colour, style=tk.ARC,width=Bond_thickness,tags=("bonds","arcs_left_slant","arcs"))
                    canvas_polygons[domain] = [arcs_left_slant[i][0], "-L-"]
        if arcs_right_slant != []:
            for i in range(len(arcs_right_slant)):
                if arcs_right_slant[i][1] == False:
                    domain = canvas.create_arc(arcs_right_slant[i][0], start=270, extent=120, fill = bond_colour, style=tk.ARC,width=Bond_thickness,tags=("bonds","arcs_right_slant","arcs"))
                    canvas_polygons[domain] = [arcs_right_slant[i][0], "-"]
                elif  arcs_right_slant[i][1] == True:
                    domain = canvas.create_arc(arcs_right_slant[i][0], start=270, extent=120, outline = linker_colour, style=tk.ARC,width=Bond_thickness,tags=("bonds","arcs_right_slant","arcs"))
                    canvas_polygons[domain] = [arcs_right_slant[i][0], "-L-"]
        for i in range(len(Linkers)):
            domain = canvas.create_line(Linkers[i], fill=linker_colour, width = Bond_thickness,tags="bonds", arrow=tk.LAST, arrowshape=Arrow_dimensions)
            canvas_polygons[domain] = [Linkers[i], "-L-"]
        for i in range(len(Hinges)):
            domain = canvas.create_line(Hinges[i][0][0], fill=hinge_colour, width = Bond_thickness,tags=("bonds","hinges"), arrow=tk.LAST, arrowshape=Arrow_dimensions)
            print(Hinges)
            if Hinges[i][1] == True:
                canvas_polygons[domain] = [Hinges[i][0][0], "-H*-"]
            else:
                canvas_polygons[domain] = [Hinges[i][0][0], "-H-"]

    #A domains
        for i in range(len(Heavy_Domains_a)):
            domain = canvas.create_polygon(Heavy_Domains_a[i], outline='#000000',fill=specificity_colours[0], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_a[i], names_Heavy_a[i]]
        for i in range(len(Light_Domains_a)):
            domain = canvas.create_polygon(Light_Domains_a[i], outline='#000000',fill=specificity_colours[1], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_a[i], names_Light_a[i]]
    #B domains
        for i in range(len(Heavy_Domains_b)):
            domain = canvas.create_polygon(Heavy_Domains_b[i], outline='#000000',fill=specificity_colours[2], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_b[i], names_Heavy_b[i]]
        for i in range(len(Light_Domains_b)):
            domain = canvas.create_polygon(Light_Domains_b[i], outline='#000000',fill=specificity_colours[3], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_b[i], names_Light_b[i]]
    #C domains
        for i in range(len(Heavy_Domains_c)):
            domain = canvas.create_polygon(Heavy_Domains_c[i], outline='#000000',fill=specificity_colours[4], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_c[i], names_Heavy_c[i]]
        for i in range(len(Light_Domains_c)):
            domain = canvas.create_polygon(Light_Domains_c[i], outline='#000000',fill=specificity_colours[5], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_c[i], names_Light_c[i]]
    #D domains
        for i in range(len(Heavy_Domains_d)):
            domain = canvas.create_polygon(Heavy_Domains_d[i], outline='#000000',fill=specificity_colours[6], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_d[i], names_Heavy_d[i]]
        for i in range(len(Light_Domains_d)):
            domain = canvas.create_polygon(Light_Domains_d[i], outline='#000000',fill=specificity_colours[7], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_d[i], names_Light_d[i]]
    #E domains
        for i in range(len(Heavy_Domains_e)):
            domain = canvas.create_polygon(Heavy_Domains_e[i], outline='#000000',fill=specificity_colours[8], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_e[i], names_Heavy_e[i]]
        for i in range(len(Light_Domains_e)):
            domain = canvas.create_polygon(Light_Domains_e[i], outline='#000000',fill=specificity_colours[9], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_e[i], names_Light_e[i]]
    #F domains
        for i in range(len(Heavy_Domains_f)):
            domain = canvas.create_polygon(Heavy_Domains_f[i], outline='#000000',fill=specificity_colours[10], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_f[i], names_Heavy_f[i]]
        for i in range(len(Light_Domains_f)):
            domain = canvas.create_polygon(Light_Domains_f[i], outline='#000000',fill=specificity_colours[11], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_f[i], names_Light_f[i]]
    #G domains
        for i in range(len(Heavy_Domains_g)):
            domain = canvas.create_polygon(Heavy_Domains_g[i], outline='#000000',fill=specificity_colours[12], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_g[i], names_Heavy_g[i]]
        for i in range(len(Light_Domains_g)):
            domain = canvas.create_polygon(Light_Domains_g[i], outline='#000000',fill=specificity_colours[13], width=1,tags="domain")
            canvas_polygons[domain] = [Light_Domains_g[i], names_Light_g[i]]
    #H domains
        for i in range(len(Heavy_Domains_h)):
            domain = canvas.create_polygon(Heavy_Domains_h[i], outline='#000000',fill=specificity_colours[14], width=1,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_h[i], names_Heavy_h[i]]
        for i in range(len(Light_Domains_h)):
            domain = canvas.create_polygon(Light_Domains_h[i], outline='#000000',fill=specificity_colours[15], width=1,tags="domain")
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
                domain = canvas.create_polygon(non_redundant_ADCs[i], outline='#000000',fill=specificity_colours[18], width=1,tags="domain")
                canvas_polygons[domain] = [non_redundant_ADCs[i],  "X"]
    #CCs
        if CCs != []:
            non_redundant_CCs = []
            non_redundant_CCs_sorted = []
            for i in range(len(CCs)):
                j = sorted(CCs[i])
                if j not in non_redundant_CCs_sorted:
                    non_redundant_CCs.append(CCs[i])
                    non_redundant_CCs_sorted.append(j)
            for i in range(len(non_redundant_CCs)):
                domain = canvas.create_polygon(non_redundant_CCs[i], outline='#000000',fill=specificity_colours[19], width=1,tags="domain")
                canvas_polygons[domain] = [non_redundant_CCs[i],  "C"]

    #Labels
        if Label_lock == True:
            for i in range(len(Label_positions)):
                x = Label_positions[i][0][0]
                y = Label_positions[i][0][1]
                Domain_Text[i] = re.sub("\_","-", Domain_Text[i])
                Domain_Text[i] = re.sub("nano|nno","",Domain_Text[i])
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
    exp = sequence_pipeline(lower_canvas)
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "fleur")
    status_label.config(text="")
    entry=textBox.get("1.0","end-1c")
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains, lower_canvas)
    render(coordinates, lower_canvas,True)


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
    lower_canvas.config(cursor = "fleur")
    status_label.config(text="")
    entry=textBox.get("1.0","end-1c")
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains, lower_canvas)
    render(coordinates, lower_canvas,True)

############################################

############################################
def save_as_png(canvas):
    fileName = "AbYdraw_export"
    file = filedialog.asksaveasfile(mode='w', filetypes=(("EPS file", "*.eps"),("PNG file", "*.png"),("JPEG file", "*.jpeg"),("All Files", "*.*") ))
    if file:
        abs_path = os.path.abspath(file.name)
        eps = canvas.postscript(file = abs_path,colormode='color',height=1000)
        file.close()
        print(abs_path)
        if "eps" not in abs_path:
            img = Image.open(str(abs_path))
            img.load(scale=10)
            if img.mode in ('P', '1'):
                img = img.convert("RGB")
            TARGET_BOUNDS = (1024, 1024)
            ratio = min(TARGET_BOUNDS[0] / img.size[0], TARGET_BOUNDS[1] / img.size[1])
            new_size = (int(img.size[0] * ratio), int(img.size[1] * ratio))
            img = img.resize(new_size, Image.ANTIALIAS)
            if "png" in abs_path:
                img.save(abs_path, "png")
            elif "jpeg" in abs_path:
                img.save(abs_path,"jpeg")



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
    global entry
    Delete_lock = False
    Bond_lock = ""
    if CLI == False:
        canvas.bind("<Button-1>", mm.select)
        canvas.bind("<B1-Motion>", mm.drag)
        canvas.bind("<ButtonRelease-1>", mm.release)
        canvas.config(cursor = "fleur")
        status_label.config(text="")
        entry=textBox.get("1.0","end-1c")
    elif CLI == True:
        entry = entry
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
        if CLI == False:
            to_save = ""
            for i in Template_File:
                to_save += (str(i)+"\n")
            f = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
            if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
                return
            f.write(to_save)
            f.close()
        elif CLI == True:
            output = open(output_name+".txt", "w+")
            for i in Template_File:
                output.write(str(i)+"\n")



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
    if CLI == False:
        String = textBox.get("1.0","end")
        Template_File.append(str("AbML String: "+String))
    elif CLI== True:
        String = entry
        Template_File.append(str("AbML String: "+String))
    save_txt_file(Template_File)




############################################
def sequence_pipeline(canvas):
    '''
    Take drawing from cavnas and generate expression to be displayed on entry box
    '''
    global all_buttons
    global Pairing_sensitivity
    for i in range(len(all_buttons)):
        all_buttons[i].config(fg="black")
    lower_canvas.config(cursor = "arrow")
    Bond_lock = ""
    Delete_lock = False
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "fleur")
    polygons_keyslist = list(canvas_polygons.keys())
    TYPE_keyslist = list(TYPE_labels.keys())
    NOTE_keyslist = list(NOTE_labels.keys())
    MOD_keyslist  = list(MOD_labels.keys())
    ANTI_keyslist = list(ANTI_labels.keys())
    LENGTH_keyslist=list(LENGTH_labels.keys())
    domains_list = lower_canvas.find_withtag("domain")
    domains_dict = {}
    for i in range(len(domains_list)):
        for j in range(len(polygons_keyslist)):
            if domains_list[i] == polygons_keyslist[j]:
                domains_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    bonds_list = lower_canvas.find_withtag("bonds")
    bonds_dict = {}
    for i in range(len(bonds_list)):
        for j in range(len(polygons_keyslist)):
            if bonds_list[i] == polygons_keyslist[j]:
                bonds_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    disulphides_list = lower_canvas.find_withtag("disulphide")
    disulphides_dict = {}
    for i in range(len(disulphides_list)):
        for j in range(len(polygons_keyslist)):
            if disulphides_list[i] == polygons_keyslist[j]:
                disulphides_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    type_list = lower_canvas.find_withtag("TYPE_labels")
    type_dict = {}
    for i in range(len(type_list)):
        for j in range(len(TYPE_keyslist)):
            if type_list[i] == TYPE_keyslist[j]:
                type_dict[j] = TYPE_labels.get(TYPE_keyslist[j])
    note_list = lower_canvas.find_withtag("NOTE_labels")
    note_dict = {}
    for i in range(len(note_list)):
        for j in range(len(NOTE_keyslist)):
            if note_list[i] == NOTE_keyslist[j]:
                note_dict[j] = NOTE_labels.get(NOTE_keyslist[j])
    mod_list = lower_canvas.find_withtag("MOD_labels")
    mod_dict = {}
    for i in range(len(mod_list)):
        for j in range(len(MOD_keyslist)):
            if mod_list[i] == MOD_keyslist[j]:
                mod_dict[j] = MOD_labels.get(MOD_keyslist[j])
    anti_list = lower_canvas.find_withtag("ANTI_labels")
    anti_dict = {}
    for i in range(len(anti_list)):
        for j in range(len(ANTI_keyslist)):
            if anti_list[i] == ANTI_keyslist[j]:
                anti_dict[j] = ANTI_labels.get(ANTI_keyslist[j])
    length_list = lower_canvas.find_withtag("LENGTH_labels")
    length_dict = {}
    for i in range(len(length_list)):
        for j in range(len(LENGTH_keyslist)):
            if length_list[i] == LENGTH_keyslist[j]:
                length_dict[j] = LENGTH_labels.get(LENGTH_keyslist[j])
    arcs_list = lower_canvas.find_withtag("arcs")
    arcs_dict = {}
    for i in range(len(arcs_list)):
        for j in range(len(polygons_keyslist)):
            if arcs_list[i] == polygons_keyslist[j]:
                arcs_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    hinges_list = lower_canvas.find_withtag("hinges")
    hinges_dict = {}
    for i in range(len(hinges_list)):
        for j in range(len(polygons_keyslist)):
            if hinges_list[i] == polygons_keyslist[j]:
                hinges_dict[j] = canvas_polygons.get(polygons_keyslist[j])

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
                    bondx2 = arc_checker(bonds_keyslist[j],bonds_dict, domains_dict)[2]
                    bondy2 = arc_checker(bonds_keyslist[j],bonds_dict, domains_dict)[3]
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
        timeout = time.time() + 2
        while continuation_found == True:
            if full_chain == [] or time.time() > timeout:
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
                            bondx1 = arc_checker(bonds_keyslist[j],bonds_dict, domains_dict)[0]
                            bondy1 = arc_checker(bonds_keyslist[j],bonds_dict, domains_dict)[1]
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
                        bondx2 = arc_checker(full_chain[-1],bonds_dict, domains_dict)[2]
                        bondy2 = arc_checker(full_chain[-1],bonds_dict, domains_dict)[3]


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
                                    bondx1 = arc_checker(bonds_keyslist[n],bonds_dict, domains_dict)[0]
                                    bondy1 = arc_checker(bonds_keyslist[n],bonds_dict, domains_dict)[1]
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
                                    bondx2 = arc_checker(full_chain[-1],bonds_dict, domains_dict)[2]
                                    bondy2 = arc_checker(full_chain[-1],bonds_dict, domains_dict)[3]
                            hingex1 = hinges_dict.get(hinges_keyslist[j])[0][0]
                            hingey1 = hinges_dict.get(hinges_keyslist[j])[0][1]

                            if -25<=(bondx2-hingex1)<=25 and -25<=(bondy2-hingey1)<=25:
                                full_chain.append(hinges_keyslist[j])
                                string.append(hinges_dict.get(hinges_keyslist[j])[1])
                        except TypeError:
                            pass




        #full_directions.append(directions)
        full_chains.append(full_chain)
        strings.append(string)

##reorder chains by specificity
    #specificities = ["a","b","c","d","e","f","g","h"]
    #categorised_chains = []
    #reordered_strings = []
    #reordered_chains = []
    #not_specified_strings = []
    #not_specified_chains  = []
    #specificities_lists_strings = [[],[],[],[],[],[],[],[]]
    #specificities_lists_chains  = [[],[],[],[],[],[],[],[]]
    #first_appearances =           [[],[],[],[],[],[],[],[]]
    #for i in range(len(strings)):
    #    for j in range(len(specificities)):
    #        if specificities[j] in str(strings[i]) and full_chains[i] not in categorised_chains:
    #            specificities_lists_strings[j].append(strings[i])
    #            specificities_lists_chains[j].append(full_chains[i])
    #            categorised_chains.append(full_chains[i])
    #            for k in range(len(strings[i])):
    #                if specificities[j] in str(strings[i][k]):
    #                    first_appearances[j].append(k)
    #                    break
    #    if full_chains[i] not in categorised_chains:
    #        not_specified_strings.append(strings[i])
    #        not_specified_chains.append(full_chains[i])
    #        categorised_chains.append(full_chains[i])

    #for i in range(len(specificities_lists_strings)):
    #    if all(x == first_appearances[i][0] for x in first_appearances[i]):
    #        specificities_lists_strings[i] = specificities_lists_strings[i]
    #        specificities_lists_chains[i] = specificities_lists_chains[i]
    #    else:
    #        specificities_lists_strings[i] = [x for _,x in sorted(zip(first_appearances[i],specificities_lists_strings[i]))]
    #        specificities_lists_chains[i] = [x for _,x in sorted(zip(first_appearances[i],specificities_lists_chains[i]))]

    #for i in range(len(specificities_lists_strings)):
    #    for j in range(len(specificities_lists_strings[i])):
    #        reordered_strings.append(specificities_lists_strings[i][j])
    #        reordered_chains.append(specificities_lists_chains[i][j])

    #for i in range(len(not_specified_strings)):
    #    reordered_strings.append(not_specified_strings[i])
    #    reordered_chains.append(not_specified_chains[i])

    #strings = reordered_strings
    #full_chains = reordered_chains

#number chains
    counter = 1
    assigned_numbers = {}
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            assigned_keyslist = list(assigned_numbers.keys())
            if str(strings[i][j]) != "-":
                if "-H-" not in strings[i][j] and "-L-" not in strings[i][j] and "-H*-" not in strings[i][j]:
                    index = full_chains[i][j]
                    coordinates = domains_dict.get(index)[0]
                    assigned_match = False
                    assigned_keyslist_check = ""
                    for k in range(len(assigned_keyslist)):
                        assigned_coordinates = assigned_numbers.get(assigned_keyslist[k])
                        if coordinates == assigned_coordinates:
                            assigned_match = True
                            assigned_keyslist_check = assigned_keyslist[k]
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
                elif "-L-" in strings[i][j]:
                    strings[i][j] = str("-L("+str(counter)+")-")
                    counter += 1
                elif "-H*-" in strings[i][j]:
                    strings[i][j] = str("-H*("+str(counter)+")-")
                    counter += 1


##Pair chains based on closeness
    paired = []
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            if ":" not in str(strings[i][j]) and "-" not in str(strings[i][j]) and "VHH" not in str(strings[i][j]):
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
                            d1x1 = min_max[0]-Pairing_sensitivity
                        else:
                            d1x2 = min_max[1]+Pairing_sensitivity
                    else:
                        d1x1 = min_max[0]-Pairing_sensitivity
                        d1x2 = min_max[1]+Pairing_sensitivity
                        d1y1 = min_max[2]
                        d1y2 = min_max[3]
                    #Search for overlapping matching domains
                    for f in range(len(domains_keyslist)):
                        if domains_keyslist[f] != index:
                            for a in range(len(full_chains)):
                                for b in range(len(full_chains[a])):
                                    try:
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
                                                        d2x1 = min_max[0]-Pairing_sensitivity
                                                    else:
                                                        d2x2 = min_max[1]+Pairing_sensitivity
                                                else:
                                                    d2x1 = min_max[0]-Pairing_sensitivity
                                                    d2x2 = min_max[1]+Pairing_sensitivity
                                                    d2y1 = min_max[2]
                                                    d2y2 = min_max[3]
                                                combinations_to_try = [[d2x1,((d2y1+d2y2)/2)],[d2x2,((d2y1+d2y2)/2)],[d2x1,d2y1],[d2x2,d2y1],[d2x2,d2y2],[d2x1,d2y2]]
                                                for g in combinations_to_try:
                                                    if d1x1 < g[0] < d1x2 and d1y1 < g[1] < d1y2:
                                                        if ("VH" in str(strings[i][j]) and "VL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("VL" in str(strings[i][j]) and "VH" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CL" in str(strings[i][j]) and "CH1" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH1" in str(strings[i][j]) and "CL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH2" in str(strings[i][j]) and "CH2" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH3" in str(strings[i][j]) and "CH3" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH4" in str(strings[i][j]) and "CH4" in str(domains_dict.get(domains_keyslist[f])[1])) or ("-H-" == str(strings[i][j]) and "-H-" == str(domains_dict.get(domains_keyslist[f])[1])) or ("-H*" in str(strings[i][j]) and "-H*-" in str(domains_dict.get(domains_keyslist[f])[1])) or ("-H*" in str(strings[i][j]) and "-H" in str(domains_dict.get(domains_keyslist[f])[1])) or ("-H" in str(strings[i][j]) and "-H*" in str(domains_dict.get(domains_keyslist[f])[1])) or ("X" in str(strings[i][j]) and "X" in str(domains_dict.get(domains_keyslist[f])[1])) or ( str(strings[i][j]) == "C" and  str(domains_dict.get(domains_keyslist[f])[1]) == "C"):
                                                            disulphide_count = 0

                                                            for y in range(len(disulphides_keyslist)):
                                                                disulphx1 = disulphides_dict.get(disulphides_keyslist[y])[0][0]
                                                                disulphy1 = disulphides_dict.get(disulphides_keyslist[y])[0][1]
                                                                disulphx2 = disulphides_dict.get(disulphides_keyslist[y])[0][2]
                                                                disulphy2 = disulphides_dict.get(disulphides_keyslist[y])[0][3]
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
                                    except IndexError:
                                        error_message = "Your drawing contains multiple antibodies that are not connected"
                                    #    raise_error(lower_canvas,error_message)

                elif ("X" in str(strings[i][j]) or str(strings[i][j])) == "C" :
                    domain_name = re.sub("\((.*?)\)","",str(strings[i][j]))
                    paired_X_domains = []
                    index = full_chains[i][j]
                    min_max = get_min_max_coordinates((domains_dict.get(index)[0]))
                    d1x1 = min_max[0]-Pairing_sensitivity
                    d1x2 = min_max[1]+Pairing_sensitivity
                    d1y1 = min_max[2]-Pairing_sensitivity
                    d1y2 = min_max[3]+Pairing_sensitivity
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
                                                        disulphx1 = disulphides_dict.get(disulphides_keyslist[y])[0][0]
                                                        disulphy1 = disulphides_dict.get(disulphides_keyslist[y])[0][1]
                                                        disulphx2 = disulphides_dict.get(disulphides_keyslist[y])[0][2]
                                                        disulphy2 = disulphides_dict.get(disulphides_keyslist[y])[0][3]

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

            elif ":" not in str(strings[i][j]) and "-H" in str(strings[i][j]):
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
                        if bonds_keyslist[f] != index and "-H" in bonds_dict.get(bonds_keyslist[f])[1]:
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
                                                        disulphx1 = disulphides_dict.get(disulphides_keyslist[y])[0][0]
                                                        disulphy1 = disulphides_dict.get(disulphides_keyslist[y])[0][1]
                                                        disulphx2 = disulphides_dict.get(disulphides_keyslist[y])[0][2]
                                                        disulphy2 = disulphides_dict.get(disulphides_keyslist[y])[0][3]
                                                        if ((d1x1 <= disulphx2 <= d1x2 and d1y1<= disulphy2 <= d1y2) and (d2x1 <= disulphx1 <= d2x2 and d2y1<= disulphy1 <= d2y2)) or ((d1x1 <= disulphx1 <= d1x2 and d1y1<= disulphy1 <= d1y2) and (d2x1 <= disulphx2 <= d2x2 and d2y1 <= disulphy2 <= d2y2)):
                                                            disulphide_count += 1
                                                    if "*" in str(strings[i][j]):
                                                        ijDomain_name = "H*"
                                                    else:
                                                        ijDomain_name = "H"
                                                    if "*" in str(strings[a][b]):
                                                        abDomain_name = "H*"
                                                    else:
                                                        abDomain_name = "H"
                                                    if disulphide_count == 0:
                                                        strings[i][j] = str("-"+ijDomain_name+"("+str(number)+":"+str(paired_number)+")-")
                                                        strings[a][b] = str("-"+abDomain_name+"("+str(paired_number)+":"+str(number)+")-")
                                                    elif disulphide_count > 0:
                                                        strings[i][j] = str("-"+ijDomain_name+"("+str(number)+":"+str(paired_number)+")"+"{"+str(disulphide_count)+"}-")
                                                        strings[a][b] = str("-"+abDomain_name+"("+str(paired_number)+":"+str(number)+")"+"{"+str(disulphide_count)+"}-")
                                                    paired.append(int(number))
                                                    paired.append(int(paired_number))
##Search for linkers connected by disulphide bonds
            elif ":" not in str(strings[i][j]) and "-L" in str(strings[i][j]):
                number =  re.findall("\((.*?)\)", str(strings[i][j]))
                number =  int(re.sub("\[|\'|\]","", str(number)))
                if number not in paired:
                    domain_name = re.sub("\((.*?)\)","",str(strings[i][j]))
                    index = full_chains[i][j]
                    d1x1 = (bonds_dict.get(index)[0][0])
                    d1x2 = (bonds_dict.get(index)[0][2])
                    if d1x1 > d1x2:
                        d1x1 = (bonds_dict.get(index)[0][2])-5
                        d1x2 = (bonds_dict.get(index)[0][0])+5
                    else:
                        d1x1 = (bonds_dict.get(index)[0][0])-5
                        d1x2 = (bonds_dict.get(index)[0][2])+5
                    d1y1 = (bonds_dict.get(index)[0][1])
                    d1y2 = (bonds_dict.get(index)[0][3])
                    for x in range(len(disulphides_keyslist)):
                        disulphx1 = disulphides_dict.get(disulphides_keyslist[x])[0][0]
                        disulphy1 = disulphides_dict.get(disulphides_keyslist[x])[0][1]
                        disulphx2 = disulphides_dict.get(disulphides_keyslist[x])[0][2]
                        disulphy2 = disulphides_dict.get(disulphides_keyslist[x])[0][3]
                        print(d1x1,disulphx1,d1x2,d1y1,disulphy1,d1y2)
                        print(d1x1,disulphx2,d1x2,d1y1,disulphy2,d1y2)
                        if (d1x1 <= disulphx2 <= d1x2 and d1y1<= disulphy2 <= d1y2) or (d1x1 <= disulphx1 <= d1x2 and d1y1<= disulphy1 <= d1y2):
                            for f in range(len(bonds_keyslist)):
                                if bonds_keyslist[f] != index:
                                    for a in range(len(full_chains)):
                                        for b in range(len(full_chains[a])):
                                            if int(full_chains[a][b]) == int(bonds_keyslist[f]) and ":" not in strings[a][b] and "L" in strings[a][b]:
                                                paired_domain = strings[a][b]
                                                paired_number =  re.findall("\((.*?)\)", str(paired_domain))
                                                paired_number =  str(re.sub("\[|\'|\]","", str(paired_number)))
                                                paired_name   = re.sub("\((.*?)\)", "", str(paired_domain))
                                                d2x1 = (bonds_dict.get(bonds_keyslist[f])[0][0])
                                                d2x2 = (bonds_dict.get(bonds_keyslist[f])[0][2])
                                                if d2x1 > d1x2:
                                                    d2x1 = (bonds_dict.get(bonds_keyslist[f])[0][2])-5
                                                    d2x2 = (bonds_dict.get(bonds_keyslist[f])[0][0])+5
                                                else:
                                                    d2x1 = (bonds_dict.get(bonds_keyslist[f])[0][0])-5
                                                    d2x2 = (bonds_dict.get(bonds_keyslist[f])[0][2])+5
                                                d2y1 = (bonds_dict.get(bonds_keyslist[f])[0][1])
                                                d2y2 = (bonds_dict.get(bonds_keyslist[f])[0][3])
                                                if ((d1x1 <= disulphx2 <= d1x2 and d1y1<= disulphy2 <= d1y2) and (d2x1 <= disulphx1 <= d2x2 and d2y1<= disulphy1 <= d2y2)) or ((d1x1 <= disulphx1 <= d1x2 and d1y1<= disulphy1 <= d1y2) and (d2x1 <= disulphx2 <= d2x2 and d2y1 <= disulphy2 <= d2y2)):
                                                    disulphide_count = 1
                                                    strings[i][j] = str("-L"+"("+str(number)+":"+str(paired_number)+")"+"{"+str(disulphide_count)+"}-")
                                                    strings[a][b] = str("-L"+"("+str(paired_number)+":"+str(number)+")"+"{"+str(disulphide_count)+"}-")
                                                    #paired.append(int(number))
                                                    #paired.append(int(paired_number))

    print(strings)
##Find comments on domains and not on domains
    comment_dicts = [type_dict, mod_dict, anti_dict, length_dict, note_dict]
    comment_lists = [type_keyslist, mod_keyslist, anti_keyslist, length_keyslist, note_keyslist]
    for i in range(len(full_chains)):
        for j in range(len(full_chains[i])):
            #print(domains_dict.get(full_chains[i][j])[0])
            domain_type = strings[i][j].split("(")[0]
            mod_domain_type = re.sub("\.|\+|\-|\@|\>","", str(domain_type))
            domain_type = re.sub("\.|\*|\+|\-|\@|\>","", str(domain_type))
            if "-" not in strings[i][j]:
                coordinates = domains_dict.get(full_chains[i][j])[0]
            elif "-" in strings[i][j]:
                coordinates = bonds_dict.get(full_chains[i][j])[0]
            min_max = get_min_max_coordinates(coordinates)
            d1x1 = min_max[0]
            d1x2 = min_max[1]
            d1y1 = min_max[2]
            d1y2 = min_max[3]
            for k in range(len(comment_lists)):
                for l in range(len(comment_lists[k])):
                    comment = comment_dicts[k].get(comment_lists[k][l])[1]
                    labelx = comment_dicts[k].get(comment_lists[k][l])[0][0]
                    labely = comment_dicts[k].get(comment_lists[k][l])[0][1]
                    #print(labelx,labely)
                    if ((d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2)) or (("X" or "C" or "-") not in mod_domain_type and mod_domain_type in str(comment) and "*" in mod_domain_type) or (("X" or "C") in mod_domain_type and mod_domain_type in str(comment)) or ((("H" or "L") in mod_domain_type and (comment[0] == "H" or comment[0] == "L")) and mod_domain_type in str(comment) and "*" in mod_domain_type):
                        note = comment_dicts[k].get(comment_lists[k][l])[1]
                        comments = note.split(",")
                        comment_to_add = ""
                        for x in range(len(comments)):
                            note1 = comments[x].split(":")[0]
                            note2 = comments[x].split(":")[1]
                            note1 = re.sub(domain_type, "",note1)
                            note1 = re.sub("^ |\* |\*","", note1)
                            noting = str("["+note1+":"+note2+"]")
                            comment_to_add +=noting
                        strings[i][j] = re.sub("\*","", strings[i][j])
                        if "MOD" in strings[i][j] or "MOD" in noting and "X" not in strings[i][j]:
                            strings[i][j] = strings[i][j].split("(")[0]+"*("+strings[i][j].split("(")[1]
                        if "-" not in strings[i][j]:
                            strings[i][j] += comment_to_add
                        elif "-" in strings[i][j]:
                            strings[i][j] = strings[i][j].split("-")[1]
                            strings[i][j] = str("-"+strings[i][j]+comment_to_add+"-")
                    else:
                        strings[i][j] = re.sub("\*","",strings[i][j])

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
    final_string = re.sub("\]\[",", ", str(final_string))
    final_string = re.sub("\-\-","-", str(final_string))
    final_string = re.sub("\-L\-\|","|",str(final_string))
    final_string = re.sub("\-\|","|",str(final_string))
    final_string = re.sub("\-$","",str(final_string))


##display string to textbox
    textBox.delete("1.0","end")
    textBox.insert("1.0",str(final_string))


def button_hover(e):
    button["bg"] = "white"
def button_hover_leave(e):
    button["bg"] = "#CFCFCF"
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

################Keypad buttons############
def update_domain_primer(domain_type,domain_charge,domain_mod,extra_mods):
    global Domain_Primer
    global Domain_Primer_Lock
    global bond_buttons
    global delete_buttons
    global comments_buttons
    global Bond_lock
    global Delete_lock
    global all_buttons
    for i in all_buttons:
        i.config(fg="black")
    Delete_lock = False
    Bond_lock = ""

    if Domain_Primer == [] and domain_type == "" and domain_charge == "" and domain_mod == "" and extra_mods == "":
        ##Turn off domain primer###
        Domain_Primer_Lock = ""
        lower_canvas.config(cursor = "fleur")
        lower_canvas.unbind("<Button-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        domain_name = ""
    elif Domain_Primer == []  and  (domain_charge != "" or  domain_mod != "" or extra_mods != "" or domain_type != ""):
        print("BOOYAA")
        ###change specificitiy###
        domain_name = str(domain_type+domain_mod+domain_charge+extra_mods)
        Domain_Primer_Lock = ""
        lower_canvas.config(cursor = "arrow")
        lower_canvas.bind("<Button-1>", mm.change_specificity)
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
    elif Domain_Primer != []:
        ###update domain_type###
        if "V" in Domain_Primer[8]:
            Domain_Primer[8] = re.sub("\.|a|b|c|d|e|f|g|h","",Domain_Primer[8])
            Domain_Primer[8] += str("."+domain_type)
        ###update KIH mods####
        Domain_Primer[5] = domain_mod
        ###update charge####
        Domain_Primer[8] = re.sub("\+|\-|\_","",Domain_Primer[8])
        if "V" in Domain_Primer[8]:
            Domain_Primer[8] = re.sub("\.",domain_charge+".",Domain_Primer[8])
        else:
            Domain_Primer[8] += domain_charge
        ##update extra mods###
        Domain_Primer[8] = re.sub("\!|\*","",Domain_Primer[8])
        if "V" in Domain_Primer[8]:
            Domain_Primer[8] = re.sub("\.",extra_mods+".",Domain_Primer[8])
        else:
            Domain_Primer[8] += extra_mods
        Domain_Primer[8] = Domain_Primer[8]
        Domain_Primer[5] = Domain_Primer[5]
        lower_canvas.bind("<Button-1>", mm.place_domain)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
        lower_canvas.config(cursor = "plus")
        domain_name = Domain_Primer[8]


        ###Update buttons###
    if  "VHH" in domain_name:
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
    if "a" in domain_type:
        a_button.config(fg="red")
    if "b" in domain_type:
        b_button.config(fg="red")
    if "c" in domain_type:
        c_button.config(fg="red")
    if "d" in domain_type:
        d_button.config(fg="red")
    if "e" in domain_type:
        e_button.config(fg="red")
    if "f" in domain_type:
        f_button.config(fg="red")
    if "g" in domain_type:
        g_button.config(fg="red")
    if "h" in domain_type:
        h_button.config(fg="red")
    if "!" in extra_mods:
        Gly_button.config(fg="red")
    if "*" in extra_mods:
        Mod_button.config(fg="red")
    if domain_mod == ">":
        KIH_knob.configure(fg="red")
    elif domain_mod == "@":
        KIH_hole.configure(fg="red")
    if domain_charge == "+":
        Positive_charge.configure(fg="red")
    elif domain_charge == "_":
        Negative_charge.configure(fg="red")

def prime_domain_button(canvas,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H,domain_name,Light,Heavy):
    global domain_type
    global domain_mod
    global extra_mods
    global domain_charge
    global Domain_Primer
    global Domain_Primer_Lock
    if Domain_Primer_Lock == "" or Domain_Primer_Lock != domain_name:
        Domain_Primer = [righthanded,slant,V,direction,X,mod,interaction,previous_H,domain_name,Light,Heavy,False,4]
        Domain_Primer_Lock = domain_name
    elif Domain_Primer_Lock == domain_name:
        Domain_Primer_Lock = ""
        Domain_Primer = []
    update_domain_primer(domain_type,domain_charge,domain_mod,extra_mods)


def domain_type_button(letter):
    global domain_type
    global domain_mod
    global extra_mods
    global domain_charge
    global Domain_Primer
    if domain_type == "" or letter not in domain_type:
        domain_type = str(letter)
    elif domain_type != "" and letter not in domain_type:
        domain_type += str(letter)
    elif letter in domain_type:
        domain_type = re.sub(letter,"",str(domain_type))
    update_domain_primer(domain_type,domain_charge,domain_mod,extra_mods)

def domain_mod_button(letter):
    global domain_type
    global domain_mod
    global extra_mods
    global domain_charge
    global Domain_Primer
    if domain_mod == "" or letter not in domain_mod:
        domain_mod = str(letter)
    elif domain_mod != "" and letter not in domain_mod:
        domain_mod += str(letter)
    elif letter in domain_mod:
        domain_mod = re.sub(letter,"",str(domain_mod))
    update_domain_primer(domain_type,domain_charge,domain_mod,extra_mods)

def extra_mod_button(letter):
    global domain_type
    global domain_mod
    global extra_mods
    global domain_charge
    global Domain_Primer
    if extra_mods == "" or letter not in extra_mods:
        extra_mods = str(letter)
    elif extra_mods != "" and letter not in extra_mods:
        extra_mods += str(letter)
    elif letter in extra_mods:
        if letter == "!":
            extra_mods = re.sub("\!","",str(extra_mods))
        elif letter == "*":
            extra_mods = re.sub("\*","",str(extra_mods))
    update_domain_primer(domain_type,domain_charge,domain_mod,extra_mods)

def domain_charge_button(letter):
    global domain_type
    global domain_mod
    global extra_mods
    global domain_charge
    global Domain_Primer
    if domain_charge == "" or letter not in domain_charge:
        domain_charge = str(letter)
    elif domain_charge != "" and letter not in domain_charge:
        domain_charge += str(letter)
    elif letter in domain_charge:
        if letter == "+":
            domain_charge = re.sub("\+","",str(domain_charge))
        elif letter == "_":
            domain_charge = re.sub("\_","",str(domain_charge))
    update_domain_primer(domain_type,domain_charge,domain_mod,extra_mods)

##########################################
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
        lower_canvas.config(cursor = "fleur")
        lower_canvas.unbind("<Button-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")


def arc_checker(key,bonds_dict, domains_dict):
    domains_keyslist = list(domains_dict.keys())
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

            min_max = get_min_max_coordinates(domains_dict.get(domains_keyslist[x])[0])
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
        lower_canvas.config(cursor = "fleur")

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
        lower_canvas.config(cursor = "fleur")


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
        lower_canvas.config(cursor = "fleur")

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
        lower_canvas.config(cursor = "fleur")

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
                if min_max[1]-min_max[0] == 70:
                    if canvas_polygons.get(polygons_keyslist[i])[0][2] > canvas_polygons.get(polygons_keyslist[i])[0][0]:
                        labelx-=5
                    elif canvas_polygons.get(polygons_keyslist[i])[0][2] < canvas_polygons.get(polygons_keyslist[i])[0][0]:
                        labelx+=5
                label  = lower_canvas.create_text([labelx,labely], text = domain_name, tags = "label")
                canvas_labels[label] = [[labelx, labely], domain_name]
                temp_label = {}
            elif domain_label == "-H-" or domain_label == "-H*-" :
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
        lower_canvas.config(cursor = "fleur")

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
        lower_canvas.config(cursor = "fleur")

def raise_error(canvas,message):
    global CLI
    if CLI == False:
        width = lower_canvas.winfo_width()
        height = lower_canvas.winfo_height()
        if width == 1 and height == 1:
            width = 1000
            height = 900
        clear_message = "click to remove error message"
        lower_canvas.delete("all")
        lower_canvas.create_text((width/2),(height/3), text = message, fill = "red")
        lower_canvas.create_text((width/2),(height/3)*2, text = clear_message)
        lower_canvas.bind("<Button-1>", mm.clear_error_message)
    elif CLI == True:
        print(message)
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
    lower_canvas.config(cursor = "fleur")
    status_label.config(text="")
    i=Library.curselection()

    index = i[0]
    entry = antibodyformats.get(formats_keyslist[index])
    textBox.insert("1.0",str(entry))
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains, lower_canvas)
    render(coordinates,lower_canvas,True)

class MouseMover():
    startcoordinates = []
    newcoordinates = []
    label = ""
    def __init__(self):
        self.item = 0; self.previous = (0, 0)
        ###Click and drag item on canvas####
    def select(self, event):
        global lower_canvas
        global canvas_labels
        global temp_label
        global canvas_polygons
        self.startcoordinates = []
        self.newcoordinates =[]
        widget =  lower_canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)

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
                    labelx_theory = (x1+x2)/2
                    labely_theory = (y1+y2)/2
                    label_text_test = re.sub("\.|\>|\@","",canvas_polygons.get(self.item)[1])
                    label_text_test = re.sub("\_","-",label_text_test)
                    print(label_text_test, label_text)
                    if (labelx_theory == labelx or labelx_theory-5 == labelx or labelx_theory+5 == labelx) and labely_theory == labely  and label_text_test==label_text:
                        del canvas_labels[label_keyslist[i]]
                        lower_canvas.delete(label_keyslist[i])
                        temp_label = {}
                        temp_label[label_keyslist[i]] = [[labelx,labely], label_text]
                return(startcoordinates)
            else:
                self.item = 0
        elif "-" in str(canvas_polygons.get(self.item)[1]):
            self.startcoordinates = [xc, yc]
            lower_canvas.config(cursor = "fleur")
            return(startcoordinates)

    def drag(self, event):
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        lower_canvas.move(self.item, xc-self.previous[0], yc-self.previous[1])
        self.previous = (xc, yc)
        coordinates = canvas_polygons.get(self.item)
        self.newcoordinates = [xc,yc]
        return(newcoordinates)
    def release(self, event):

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

            if self.newcoordinates != self.startcoordinates[0]:
                arcs_list = lower_canvas.find_withtag("arcs")
                min_max = get_min_max_coordinates(coordinates)
                x1 = min_max[0]
                x2 = min_max[1]
                y1 = min_max[2]
                y2 = min_max[3]

                for i in range(len(canvas_keyslist)):
                    if "-" in str(canvas_polygons.get(canvas_keyslist[i])[1]):
                        domain_coordinates = (canvas_polygons.get(canvas_keyslist[i])[0])
                        bondx1 = domain_coordinates[0]
                        bondy1 = domain_coordinates[1]
                        bondx2 = domain_coordinates[2]
                        bondy2 = domain_coordinates[3]
                        if canvas_keyslist[i] in arcs_list:
                            domains_list = lower_canvas.find_withtag("domain")
                            domains_dict = {}
                            for x in range(len(domains_list)):
                                for y in range(len(canvas_keyslist)):
                                    if domains_list[x] == canvas_keyslist[y]:
                                        domains_dict[y] = canvas_polygons.get(canvas_keyslist[y])
                            bonds_list = lower_canvas.find_withtag("bonds")
                            bonds_dict = {}
                            for x in range(len(bonds_list)):
                                for y in range(len(canvas_keyslist)):
                                    if bonds_list[x] == canvas_keyslist[y]:
                                        bonds_dict[y] = canvas_polygons.get(canvas_keyslist[y])
                            bondx1 = arc_checker(canvas_keyslist[i],bonds_dict, domains_dict)[0]
                            bondy1 = arc_checker(canvas_keyslist[i],bonds_dict, domains_dict)[1]
                            bondx2 = arc_checker(canvas_keyslist[i],bonds_dict, domains_dict)[2]
                            bondy2 = arc_checker(canvas_keyslist[i],bonds_dict, domains_dict)[3]
                        if x1< bondx2 <x2 and y1 < bondy2  < y2 or  x1< bondx1 <x2 and y1 < bondy1  < y2:
                            lower_canvas.delete(canvas_keyslist[i])
                            del canvas_polygons[canvas_keyslist[i]]
            if Label_lock == True:
                temp_label_key = list(temp_label.keys())
                if len(temp_label_key) >0:
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
                if x2-x1 == 70:
                    if new_coordinates[2] > new_coordinates[0]:
                        labelx2-=5
                    elif new_coordinates[2] < new_coordinates[0]:
                        labelx2+=5
                label_name = temp_label.get(temp_label_key[0])[1]
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
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        polygons_keyslist = list(canvas_polygons.keys())
        labels_keyslist = list(canvas_labels.keys())
        TYPE_keyslist = list(TYPE_labels.keys())
        NOTE_keyslist = list(NOTE_labels.keys())
        MOD_keyslist = list(MOD_labels.keys())
        ANTI_keyslist = list(ANTI_labels.keys())
        LENGTH_keyslist = list(LENGTH_labels.keys())
        domain_self_item = False
        for i in range(len(polygons_keyslist)):
            domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            if x1< xc <x2 and y1 < yc < y2:
                self.item = polygons_keyslist[i]
                domain_self_item = True

        if domain_self_item == False:
            for i in range(len(labels_keyslist)):
                label_coordinates = (canvas_labels.get(labels_keyslist[i])[0])
                lab_x1 = min_max[0]
                lab_y1 = min_max[1]
                if lab_x1 == xc and lab_y1 == yc:
                    for j in range(len(polygons_keyslist)):
                        domain_coordinates = (canvas_polygons.get(polygons_keyslist[j])[0])
                        min_max = get_min_max_coordinates(domain_coordinates)
                        x1 = min_max[0]
                        x2 = min_max[1]
                        y1 = min_max[2]
                        y2 = min_max[3]
                        if x1 <= lab_x1 <= x2 and y1 <= lab_y1 <= y2:
                            self.item = polygons_keyslist[j]
                            domain_self_item = True
        if domain_self_item == False:

            self.item = widget.find_closest(xc, yc)[0]
            if self.item not in polygons_keyslist and self.item not in labels_keyslist:
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
            elif self.item in polygons_keyslist and "disulphide" in canvas_polygons.get(self.item)[1]:
                lower_canvas.delete(self.item)
                del canvas_polygons[self.item]
            else:
                self.item = None


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
        widget = lower_canvas                       # Get handle to canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        global canvas_polygons
        global domain_type
        global domain_mod
        global extra_mods
        global domain_charge
        global canvas_labels
        global temp_label
        global specificity_colours
        domain_self_item = False
        label_keyslist = list(canvas_labels.keys())
        polygons_keyslist = list(canvas_polygons.keys())
        for i in range(len(polygons_keyslist)):
            domain_coordinates = (canvas_polygons.get(polygons_keyslist[i])[0])
            min_max = get_min_max_coordinates(domain_coordinates)
            x1 = min_max[0]
            x2 = min_max[1]
            y1 = min_max[2]
            y2 = min_max[3]
            if x1< xc <x2 and y1 < yc < y2:
                self.item = polygons_keyslist[i]
                domain_self_item = True
        if domain_self_item == True:
            domain_name = canvas_polygons.get(self.item)[1]
            clean_domain_name = re.sub("\.a|\.b|\.c|\.d|\.e|\.f|\.g|\.h|a|b|c|d|e|f|g|h","",domain_name)
            clean_domain_name = re.sub("\@|\>|\+|\-|\_|\!|\*| ","",clean_domain_name)
            if domain_type != "":
                domain_type_to_add = re.sub("\.","",domain_type)
                domain_type_to_add = str("."+domain_type_to_add)
            elif domain_type == "":
                if "." not in domain_name:
                    domain_type_to_add = ""
                    for s in range(len(domain_name)):
                        if domain_name[s].islower():
                            domain_type_to_add += domain_name[s]
                else:
                    domain_type_to_add = domain_name.split(".")[1]

            if "+" in str(domain_name) or "_" in str(domain_name):
                if "+" in domain_charge and str(domain_name):
                    domain_charge_to_add = re.sub("\+","",domain_charge)
                elif "_" in domain_charge and str(domain_name):
                    domain_charge_to_add = re.sub("\_","",domain_charge)
            else:
                domain_charge_to_add = domain_charge

            if "@" in str(domain_name) or ">" in str(domain_name):
                domain_mod_to_add = re.sub(domain_mod,"",domain_mod)
            else:
                domain_mod_to_add = domain_mod

            if "*" in str(domain_name) or "!" in str(domain_name):
                if "*" in extra_mods:
                    extra_mods_to_add = re.sub("\*","",extra_mods)
                if "!" in extra_mods:
                    extra_mods_to_add = re.sub("\!","",extra_mods)
            else:
                extra_mods_to_add = extra_mods

            new_domain_name= str(clean_domain_name+domain_mod_to_add+extra_mods_to_add+domain_charge_to_add+domain_type_to_add)
            new_domain_name = re.sub(" ","",new_domain_name)
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
                    new_display_name = re.sub("\_","-", new_domain_name)
                    new_display_name = re.sub("\.| ","", new_display_name)
                    label  = lower_canvas.create_text([labelx,labely], text = new_display_name, tags = "label")
                    canvas_labels[label] = [[labelx,labely], new_display_name]
                    temp_label[label] = [[labelx,labely], new_domain_name]

        elif domain_self_item == False:
            return
    def change_orientation(self,event):
        self.startcoordinates = []
        self.newcoordinates =[]
        widget = lower_canvas                       # Get handle to canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
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
        x1y1 = (canvas_polygons.get(self.item)[0][2])
        x1y2 = (canvas_polygons.get(self.item)[0][4])
        if x1y1 == x1y2 and x1< xc <x2 and y1 < yc < y2:
            new_coordinates = []
            if "X" not in domain_name:
                for i in range(len(coordinates)):
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
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=heavy_colour, width=1, tags="domain")
            elif "VL" in domain_name or "CL" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=light_colour, width=1, tags="domain")
            elif "X" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=X_colour, width=1, tags="domain")
            elif "C" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=C_colour, width=1, tags="domain")


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


        else:
            self.item = 0
    ###Place Domains ####
    def place_domain(self,event):
        global Domain_Primer
        global Label_lock
        global specificity_colours
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        startx = xc
        starty = yc-40
        Domain_Primer[5] = re.sub("\+|\-|\_","",Domain_Primer[5])
        if Domain_Primer[5] == "@" or Domain_Primer[5] == ">":
            global domain_direction
            Domain_Primer[3] = "innie"
        All_positions_and_chains = {}
        domaincoordinates = domainmaker(All_positions_and_chains,startx,starty,Domain_Primer[0],Domain_Primer[1],Domain_Primer[2],Domain_Primer[3],Domain_Primer[4],Domain_Primer[5],Domain_Primer[6],Domain_Primer[7],Domain_Primer[11])
        domain_name = re.sub("nano|nno","",Domain_Primer[8])
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
            domain = lower_canvas.create_polygon(domaincoordinates[0], outline='#000000',fill=heavy_colour, width=1, tags="domain")
        elif Domain_Primer[9] == True:
            domain = lower_canvas.create_polygon(domaincoordinates[0], outline='#000000',fill=light_colour, width=1, tags="domain")
        canvas_polygons[domain] = [domaincoordinates[0], domain_name]
        if Label_lock == True:
            domain_name = re.sub("\.|@|>","",domain_name)
            domain_name = re.sub("\_","-",domain_name)
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
        global Type_options
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "TYPE_labels")
            TYPE_labels[label] = [[xc,yc], entry]
        else:
            if Type_list_header.get() != "TYPE":
                label = lower_canvas.create_text(xc,yc, text = Type_list_header.get(), tags = "MOD_labels")
                MOD_labels[label] = [[xc,yc], Type_list_header.get()]
    def place_note_label(self,event):
        global NOTE_labels
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        entry=str("NOTE:"+CustomLabelEntry.get("1.0","end-1c"))
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "NOTE_labels")
            NOTE_labels[label] = [[xc,yc], entry]

    def place_mod_label(self,event):
        global MOD_labels
        global Mod_options
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "MOD_labels")
            MOD_labels[label] = [[xc,yc], entry]
        else:
            if mod_list_header.get() != "MOD":
                label = lower_canvas.create_text(xc,yc, text = mod_list_header.get(), tags = "MOD_labels")
                MOD_labels[label] = [[xc,yc], mod_list_header.get()]

    def place_anti_label(self,event):
        global ANTI_labels
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        entry=str("ANTI:"+CustomLabelEntry.get("1.0","end-1c"))
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "ANTI_labels")
            ANTI_labels[label] = [[xc,yc], entry]
    def place_length_label(self,event):
        global LENGTH_labels
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        entry=str("LENGTH:"+CustomLabelEntry.get("1.0","end-1c"))
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "LENGTH_labels")
            LENGTH_labels[label] = [[xc,yc], entry]
    ###Draw and drag bonds###
    def start_bond(self,event):
        lower_canvas.delete("draggable_line")
        # Convert screen coordinates to canvas coordinates
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        self.startcoordinates = [xc, yc]
        return(self.startcoordinates)
    def drag_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#000000', width=Bond_thickness,tags = "draggable_line", arrow=tk.LAST , arrowshape=Bond_Arrows)
        self.newcoordinates = [x2,y2]
        return(newcoordinates)
    def release_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[1], width=Bond_thickness, tags=("bonds","connector"), arrow=tk.LAST , arrowshape=Bond_Arrows)
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates   = []
    def release_Hinge_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-H-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[2], width=Bond_thickness, tags=("bonds","hinge"), arrow=tk.LAST , arrowshape=Bond_Arrows)
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordiantes   = []
    def release_Linker_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-L-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[3], width=Bond_thickness, tags=("bonds","linker"), arrow=tk.LAST , arrowshape=Bond_Arrows)
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates =[]
    def drag_disulphide_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[0], width=Bond_thickness,tags = "draggable_line", arrow=tk.BOTH, arrowshape=Bond_Arrows)
        self.newcoordinates = [x2,y2]
    def drag_Hinge_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[2], width=Bond_thickness,tags = "draggable_line", arrow=tk.LAST , arrowshape=Bond_Arrows)
        self.newcoordinates = [x2,y2]
    def drag_Linker_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        widget = lower_canvas
        xc = widget.canvasx(event.x); yc = widget.canvasy(event.y)
        x1 = self.startcoordinates[0]
        y1 = self.startcoordinates[1]
        x2 = xc
        y2 = yc
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[3], width=Bond_thickness,tags = "draggable_line", arrow=tk.LAST , arrowshape=Bond_Arrows)
        self.newcoordinates = [x2,y2]

    def release_Disulphide_bond(self,event):
        global Bond_Arrows
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-disulphide-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill=bond_colours[0], width=Bond_thickness, tags="disulphide", arrow=tk.BOTH, arrowshape=Bond_Arrows)
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates =[]

    def clear_error_message(self,event):
        global canvas_polygons
        global canvas_labels
        global temp_label
        global TYPE_labels
        global NOTE_labels
        global MOD_labels
        global ANTI_labels
        global LENGTH_labels
        global specificity_colours
        lower_canvas.delete("all")
        polygons_keyslist = list(canvas_polygons.keys())
        for i in range(len(polygons_keyslist)):
            domain_coordinates = canvas_polygons.get(polygons_keyslist[i])[0]
            domain_name = canvas_polygons.get(polygons_keyslist[i])[1]
            if "-" not in str(domain_name):
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
                    domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=heavy_colour, width=1, tags="domain")
                else:
                    domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=light_colour, width=1, tags="domain")


                if Label_lock == True:
                    domain_name = re.sub("\.|@|>","",domain_name)
                    min_max = get_min_max_coordinates(domain_coordinates)
                    labelx = get_min_max_coordinates(domain_coordinates)[4]
                    labely = get_min_max_coordinates(domain_coordinates)[5]
                    if min_max[1]-min_max[0] == 70:
                        if domain_coordinates[2] > domain_coordinates[0]:
                            labelx-=5
                        elif domain_coordinates[2] < domain_coordinates[0]:
                            labelx+=5
                    label  = lower_canvas.create_text(labelx,labely, text = str(domain_name), tags = "label")
                    canvas_labels[label] = [[labelx,labely], domain_name]
            elif domain_name == "-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=bond_colour, width=Bond_thickness, tags="bonds")
            elif domain_name == "-H-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=hinge_colour, width=Bond_thickness, tags="bonds")
            elif domain_name == "-L-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=linker_colour, width=Bond_thickness, tags="bonds")
            elif domain_name == "-disulphide-":####
                domain = lower_canvas.create_line(domain_coordinates, fill=disulphide_colour, width=Bond_thickness, tags=("disulphide","bonds"))
        for i in range(len(list(TYPE_labels.keys()))):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "TYPE_labels")
        for i in range(len(list(NOTE_labels.keys()))):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "NOTE_labels")
        for i in range(len(list(MOD_labels.keys()))):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "MOD_labels")
        for i in range(len(list(ANTI_labels.keys()))):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "ANTI_labels")
        for i in range(len(list(LENGTH_labels.keys()))):
            label = lower_canvas.create_text(domain_coordinates, text = entry, tags = "LENGTH_labels")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        lower_canvas.config(cursor = "fleur")



################Domain Drawer######################
def domainmaker(All_positions_and_chains,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H, Build_up):
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



    #if V == True and direction == "constant":
    #    if righthanded == False:

    global Show_Leucine_Zippers
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

    elif (X == True and mod != "Leucine") or (X == True and mod == "Leucine" and Show_Leucine_Zippers == False):
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
    elif mod == "Leucine" and Show_Leucine_Zippers == True:
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

    if (slant == True or mod == "Leucine" or previous_H == True):# and X == False:
        top_bond    = [firstx,firsty+2]
    elif Build_up == True or X == True:
        top_bond    = [firstx,firsty+20]
    else:
        top_bond    = [firstx,firsty+2]

    if mod=="Leucine" and Show_Leucine_Zippers == True:
        bottom_bond = [ninthx,ninthy]
    else:
        bottom_bond = [fourthx, fourthy-2]
    try:
        x12 = get_min_max_coordinates(coordinates)[4]
        y12 = get_min_max_coordinates(coordinates)[5]
    except ValueError:
        x12 = (firstx)
        y12 = (firsty+thirdy)/2
    if slant == True and righthanded == False and (mod != "Leucine" or Show_Leucine_Zippers == False):
        Labelbond   = [firstx+20, y12]
    elif slant == True and righthanded == True and (mod != "Leucine" or Show_Leucine_Zippers == False):
        Labelbond   = [firstx-20, y12]
    elif mod== "Leucine" and righthanded == False and Show_Leucine_Zippers == True:
        Labelbond   = [firstx-30, (firsty+ninthy)/2]
    elif mod== "Leucine" and righthanded == True and Show_Leucine_Zippers == True:
        Labelbond   = [firstx+30, (firsty+ninthy)/2]
    else:
        Labelbond   = [x12, y12]
    return(coordinates, bottom_bond,top_bond,Labelbond)



###############Main programme#######################

root = tk.Tk()
if len(sys.argv) < 2:
    CLI = False
    HEIGHT = root.winfo_screenheight()
    WIDTH  = root.winfo_screenwidth()
    def browseFiles():
        global textBox
        filename = filedialog.askopenfilename(initialdir = "/",title = "Open File",filetypes = (("Text files","*.txt"),("all files","*.*")))
        entry = Get_input(filename)
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
        global all_buttons
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
                    domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=heavy_colour, width=1, tags="domain")
                elif "L" in domain_name:
                    domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=light_colour, width=1, tags="domain")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
                if Label_lock == True:
                    domain_name = re.sub("\.|@|>","",domain_name)
                    labelx = get_min_max_coordinates(domain_coordinates)[4]
                    labely = get_min_max_coordinates(domain_coordinates)[5]
                    if get_min_max_coordinates(domain_coordinates)[1]-get_min_max_coordinates(domain_coordinates)[0] == 70 :
                        if domain_coordinates[2] > domain_coordinates[0]:
                            labelx-=5
                        elif domain_coordinates[2] < domain_coordinates[0]:
                            labelx+=5
                    label  = lower_canvas.create_text(labelx,labely, text = str(domain_name), tags = "label")
                    canvas_labels[label] = [[labelx,labely], domain_name]
            elif domain_name == "-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=bond_colour, width=Bond_thickness, tags="bonds")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            elif domain_name == "-H-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=hinge_colour, width=Bond_thickness, tags="bonds")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            elif domain_name == "-L-" :####
                domain = lower_canvas.create_line(domain_coordinates, fill=linker_colour, width=Bond_thickness, tags="bonds")
                canvas_polygons[domain] = [domain_coordinates, domain_name]
            elif domain_name == "-disulphide-":####
                domain = lower_canvas.create_line(domain_coordinates, fill=disulphide_colour, width=Bond_thickness, tags=("disulphide","bonds"))
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
                        domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=heavy_colour, width=1, tags="domain")
                    elif "L" in domain_name:
                        domain = lower_canvas.create_polygon(domain_coordinates, outline='#000000',fill=light_colour, width=1, tags="domain")
                    canvas_polygons[domain] = [domain_coordinates, domain_name]
                    if Label_lock == True:
                        domain_name = re.sub("\.|@|>","",domain_name)
                        labelx = get_min_max_coordinates(domain_coordinates)[4]
                        labely = get_min_max_coordinates(domain_coordinates)[5]
                        if get_min_max_coordinates(domain_coordinates)[1]-get_min_max_coordinates(domain_coordinates)[0] == 70:
                            if domain_coordinates[2] > domain_coordinates[0]:
                                labelx-=5
                            elif domain_coordinates[2] < domain_coordinates[0]:
                                labelx+=5
                        label  = lower_canvas.create_text(labelx,labely, text = str(domain_name), tags = "label")
                        canvas_labels[label] = [[labelx,labely], domain_name]
                elif domain_name == "-" :####
                    domain = lower_canvas.create_line(domain_coordinates, fill=bond_colour, width=Bond_thickness, tags="bonds")
                    canvas_polygons[domain] = [domain_coordinates, domain_name]
                elif domain_name == "-H-" :####
                    domain = lower_canvas.create_line(domain_coordinates, fill=hinge_colour, width=Bond_thickness, tags="bonds")
                    canvas_polygons[domain] = [domain_coordinates, domain_name]
                elif domain_name == "-L-" :####
                    domain = lower_canvas.create_line(domain_coordinates, fill=linker_colour, width=Bond_thickness, tags="bonds")
                    canvas_polygons[domain] = [domain_coordinates, domain_name]
                elif domain_name == "-disulphide-":####
                    domain = lower_canvas.create_line(domain_coordinates, fill=disulphide_colour, width=Bond_thickness, tags=("disulphide","bonds"))
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
    textBox = tk.Text(frame, font=20)
    textBox.place(relx=0.01,rely = 0.65, relwidth=0.4,relheight=0.15)
    ###Option box
    frame2 = tk.Frame(frame, bg = '#D3D3D3')
    frame2.place(relx=0.01, rely = 0.05,relheight = 0.35, relwidth = 0.4)
    ##SEQUENCE PIPELINE VARIABLES
    Domain_Primer = []
    Domain_Primer_Lock = ""
    domain_type = ""
    domain_mod  = ""
    extra_mods  = ""
    domain_charge=""
    Pairing_sensitivity = 30
    Bond_thickness = 2
    domain_direction = "constant"
    H_Labels = False
    L_Labels = False
    Show_Leucine_Zippers = False
    Bond_Arrows=(0,0,0)
    ######################BUTTONS
    a_button= tk.Button(frame2,text="a",bg = "#CFCFCF", command =lambda: domain_type_button("a"))
    a_button.place(relx = 0.61, rely = 0.21, relheight = 0.1, relwidth=0.05)
    b_button= tk.Button(frame2,text="b",bg = "#CFCFCF", command =lambda: domain_type_button("b"))
    b_button.place(relx = 0.66, rely = 0.21, relheight = 0.1, relwidth=0.05)
    c_button= tk.Button(frame2,text="c",bg = "#CFCFCF", command =lambda: domain_type_button("c"))
    c_button.place(relx = 0.71, rely = 0.21, relheight = 0.1, relwidth=0.05)
    d_button= tk.Button(frame2,text="d",bg = "#CFCFCF", command =lambda: domain_type_button("d"))
    d_button.place(relx = 0.76, rely = 0.21, relheight = 0.1, relwidth=0.05)
    e_button= tk.Button(frame2,text="e",bg = "#CFCFCF", command =lambda: domain_type_button("e"))
    e_button.place(relx = 0.61, rely = 0.31, relheight = 0.1, relwidth=0.05)
    f_button= tk.Button(frame2,text="f",bg = "#CFCFCF", command =lambda: domain_type_button("f"))
    f_button.place(relx = 0.66, rely = 0.31, relheight = 0.1, relwidth=0.05)
    g_button= tk.Button(frame2,text="g",bg = "#CFCFCF", command =lambda: domain_type_button("g"))
    g_button.place(relx = 0.71, rely = 0.31, relheight = 0.1, relwidth=0.05)
    h_button= tk.Button(frame2,text="h",bg = "#CFCFCF", command =lambda: domain_type_button("h"))
    h_button.place(relx = 0.76, rely = 0.31, relheight = 0.1, relwidth=0.05)
    Mod_button = tk.Button(frame2,text="*",bg = "#CFCFCF", command =lambda: extra_mod_button("*"))
    Mod_button.place(relx = 0.61, rely = 0.41, relheight = 0.1, relwidth=0.05)
    Gly_button = tk.Button(frame2,text="!",bg = "#CFCFCF", command =lambda: extra_mod_button("!"))
    Gly_button.place(relx = 0.66, rely = 0.41, relheight = 0.1, relwidth=0.05)
    KIH_knob= tk.Button(frame2,text=">",bg = "#CFCFCF", command =lambda: domain_mod_button(">"))
    KIH_knob.place(relx = 0.71, rely = 0.41, relheight = 0.1, relwidth=0.05)
    KIH_hole= tk.Button(frame2,text="@",bg = "#CFCFCF", command =lambda: domain_mod_button("@"))
    KIH_hole.place(relx = 0.76, rely = 0.41, relheight = 0.1, relwidth=0.05)
    Positive_charge= tk.Button(frame2,text="+",bg = "#CFCFCF", command =lambda: domain_charge_button("+"))
    Positive_charge.place(relx = 0.61, rely = 0.51, relheight = 0.1, relwidth=0.05)
    Negative_charge= tk.Button(frame2,text="-",bg = "#CFCFCF", command =lambda: domain_charge_button("_"))
    Negative_charge.place(relx = 0.66, rely = 0.51, relheight = 0.1, relwidth=0.05)
    LengthLabelButton = tk.Button(frame2,text="LENGTH", bg="#CFCFCF", command=lambda: SelectCommentTypeButton("LENGTH"))
    LengthLabelButton.place(relx = 0.71,rely = 0.51, relheight = 0.1, relwidth= 0.1)
    NoteLabelButton = tk.Button(frame2,text="NOTE", bg="#CFCFCF", command=lambda: SelectCommentTypeButton("NOTE"))
    NoteLabelButton.place(relx = 0.61,rely = 0.61, relheight = 0.1, relwidth= 0.1)
    TypeLabelButton = tk.Button(frame2,text="TYPE", bg="#CFCFCF", command=lambda: SelectCommentTypeButton("TYPE"))
    TypeLabelButton.place(relx = 0.71,rely = 0.61, relheight = 0.1, relwidth= 0.1)
    AntiLabelButton = tk.Button(frame2,text="ANTI", bg="#CFCFCF", command=lambda: SelectCommentTypeButton("ANTI"))
    AntiLabelButton.place(relx = 0.61,rely = 0.71, relheight = 0.1, relwidth= 0.1)
    ModLabelButton = tk.Button(frame2,text="MOD", bg="#CFCFCF", command=lambda: SelectCommentTypeButton("MOD"))
    ModLabelButton.place(relx = 0.71,rely = 0.71, relheight = 0.1, relwidth= 0.1)
    nanobody_button = tk.Button(frame2,text="VHH",bg = "#CFCFCF", command =lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"Single_Fv_Chain",False,domain_mod,"","",str("VHH"+domain_mod+"."+domain_type),False,True))
    nanobody_button.place(relx = 0.61, rely = 0.01, relheight = 0.2, relwidth=0.2)
    ###Insert bonds buttons ###
    ##Col1
    InsertVHDomainButton= tk.Button(frame2,text="VH",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"innie",False,domain_mod,"","",str("VH"),False,True))
    InsertVHDomainButton.place(relx = 0.21, rely = 0.01, relheight = 0.2, relwidth=0.2)
    InsertCH1DomainButton= tk.Button(frame2,text="CH1",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH1"),False,True))
    InsertCH1DomainButton.place(relx = 0.21, rely = 0.21, relheight = 0.2, relwidth=0.2)
    InsertCH2DomainButton= tk.Button(frame2,text="CH2",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH2"),False,True))
    InsertCH2DomainButton.place(relx = 0.21, rely = 0.41, relheight = 0.2, relwidth=0.2)
    InsertCH3DomainButton= tk.Button(frame2,text="CH3",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH3"),False,True))
    InsertCH3DomainButton.place(relx = 0.21, rely = 0.61, relheight = 0.2, relwidth=0.2)
    ##Col2
    InsertVLDomainButton= tk.Button(frame2,text="VL",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"outie",False,domain_mod,"","",str("VL"),True,False))
    InsertVLDomainButton.place(relx = 0.41, rely = 0.01, relheight = 0.2, relwidth=0.2)
    InsertCLDomainButton= tk.Button(frame2,text="CL",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CL"),True,False))
    InsertCLDomainButton.place(relx = 0.41, rely = 0.21, relheight = 0.2, relwidth=0.2)
    InsertCH4DomainButton= tk.Button(frame2,text="CH4",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH4"),False,True))
    InsertCH4DomainButton.place(relx = 0.41, rely = 0.41, relheight = 0.2, relwidth=0.2)
    InsertXDomainButton= tk.Button(frame2,text="X",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,True,domain_mod,"","",str("X"),True,False))
    InsertXDomainButton.place(relx = 0.41, rely = 0.61, relheight = 0.2, relwidth=0.1)
    InsertCDomainButton= tk.Button(frame2,text="C",bg = "#CFCFCF", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,True,domain_mod,"","",str("C"),True,False))
    InsertCDomainButton.place(relx = 0.51, rely = 0.61, relheight = 0.2, relwidth=0.1)

    ###Drag and pull bonds###
    Bond_lock = ""
    InsertBondButton= tk.Button(frame2,text="Connect",bg = "#CFCFCF", command=lambda: bond_drag_button(lower_canvas,"-","bond"))
    InsertBondButton.place(relx = 0.01, rely = 0.01, relheight = 0.2, relwidth=0.2)
    InsertDBondButton= tk.Button(frame2,text="Disulphide",bg = "#CFCFCF", command=lambda: disulphide_bond_button(lower_canvas,"-","disulphide"))
    InsertDBondButton.place(relx = 0.01, rely = 0.21, relheight = 0.2, relwidth=0.2)
    InsertLHingeButton= tk.Button(frame2,text="Hinge",bg = "#CFCFCF", command=lambda: Hinge_bond_button(lower_canvas,"-H-","hinge"))
    InsertLHingeButton.place(relx = 0.01, rely = 0.41, relheight = 0.2, relwidth=0.2)
    InsertLinkerButton= tk.Button(frame2,text="Linker",bg = "#CFCFCF", command=lambda: Linker_bond_button(lower_canvas,"-L-","linker"))
    InsertLinkerButton.place(relx = 0.01, rely = 0.61, relheight = 0.2, relwidth=0.2)

    ###Delete/clear###
    Delete_lock = False
    Label_lock = True
    CustomLabelLock = ""
    InsertDelAllButton = tk.Button(frame2,text="Clear All",bg="#CFCFCF", command=lambda: delete_all_button(lower_canvas))
    InsertDelAllButton.place(relx = 0.81,rely = 0.01, relheight = 0.2, relwidth= 0.18)
    InsertDelClickButton = tk.Button(frame2,text="Delete",bg="#CFCFCF", command=lambda: delete_button(lower_canvas))
    InsertDelClickButton.place(relx = 0.81,rely = 0.21, relheight = 0.2, relwidth= 0.18)
    Labels_buttons = tk.Button(frame2,text="Labels", bg="#CFCFCF", command=lambda: labels_button(lower_canvas))
    Labels_buttons.place(relx = 0.81,rely = 0.41, relheight = 0.2, relwidth= 0.18)
    Type_options = ["TYPE","TYPE:ZIPPER", "TYPE:FUSION", "TYPE:OTHER"]
    Type_list_header = tk.StringVar(frame2)
    Type_list_header.set(Type_options[0])
    Type_menu = tk.OptionMenu(frame2, Type_list_header, *Type_options)
    Type_menu.place(relx=0.01,rely = 0.83, relwidth=0.33,relheight=0.15)
    Mod_options = ["MOD", "MOD:ENHANCEFCRN","MOD:ENHANCEADCC",    "MOD:STRANDEXCHANGE","MOD:DISULPHIDE",    "MOD:DISULFIDE",    "MOD:PI",    "MOD:CONJUGATION",   "MOD:HEXAMER",    "MOD:NOFCGR",    "MOD:NOPROTEINA","MOD:NOOX",    "MOD:NOADCC",    "MOD:NOCDC",    "MOD:NOADCP",    "MOD:NOADCCCDC",    "MOD:NOGLYCOS","MOD:NOADE",    "MOD:NOAGG",    "MOD:NOPROT",    "MOD:REMCYS",    "MOD:STABILIZATION",    "MOD:AFFINITY",    "MOD:OTHER"]
    mod_list_header = tk.StringVar(frame2)
    mod_list_header.set(Mod_options[0])
    Mod_menu = tk.OptionMenu(frame2, mod_list_header, *Mod_options)
    Mod_menu.place(relx=0.34,rely = 0.83, relwidth=0.33,relheight=0.15)



    CustomLabelButton = tk.Button(frame2,text="Comment", bg="#CFCFCF", command=lambda: CommentLabelButton_function(lower_canvas))
    CustomLabelButton.place(relx = 0.81,rely = 0.61, relheight = 0.2, relwidth= 0.18)
    CustomLabelEntry = tk.Text(frame2, font=20)
    CustomLabelEntry.place(relx=0.67,rely = 0.83, relwidth=0.31,relheight=0.15)
    ##Status bar###
    status_label = tk.Label(root, text='', bd=1)
    status_label.place(rely = 0.98, relheight = 0.02, relwidth = 1)

    ##Big button
    button = tk.Button(frame, text = "Get Structure", bg = "#CFCFCF", font=20, command=lambda: render_pipeline(lower_canvas))
    button.place(relx=0.05,rely=0.825,relheight=0.1, relwidth=0.15)
    button.bind("<Enter>", button_hover)
    button.bind("<Leave>", button_hover_leave)
    button = tk.Button(frame, text = "Get AbML", bg = "#CFCFCF", font=20, command=lambda: sequence_pipeline(lower_canvas))
    button.place(relx=0.2,rely=0.825,relheight=0.1, relwidth=0.15)
    button.bind("<Enter>", button_hover)
    button.bind("<Leave>", button_hover_leave)
    button = tk.Button(frame, text = "Tidy", bg = "#CFCFCF", font=20, command=lambda: sequence_render_pipeline(lower_canvas))
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
    "Knobs in Holes":"VH.a(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3*@(5:12) | VL.a(6:1)-CL(7:2){1} | VH.b(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3*>(12:5) | VL.b(13:8)-CL(14:9){1}",
    "Orthagonal Fab":"VH*+.a(1:6)-CH1*>(2:7){1}-H(3:10){2}-CH2(4:11)-CH3*@(5:12)|VL_.a(6:1)-CL*@(7:2){1}|VH.b(8:13)-CH1*>(9:14){1}-H(10:3){2}-CH2(11:4)-CH3*>(12:5)|VL.b(13:8)-CL*@(14:9){1}",
    "IgG-H-scFV":"VH.a(1:19)-CH1(2:20){1}-H(3:12){2}-CH2(4:13)-CH3(5:14)-L(6)-VH.b(7:9)-L(8)-VL.b(9:7)|VH.a(10:21)-CH1(11:22){1}-H(12:3){2}-CH2(13:4)-CH3(14:5)-L(15)-VH.b(16:18)-L(17)-VL.b(18:16)|VL.a(19:1)-CL(20:2){1}|VL.a(21:10)-CL(22:11){1}",
    "IgG-L-scFV":"VH.a(1:11)-CH1(2:12){1}-H(3:8){2}-CH2(4:9)-CH3(5:10)|VH.a(6:17)-CH1(7:18){1}-H(8:3){2}-CH2(9:4)-CH3(10:5)|VL.a(11:1)-CL(12:2){1}-L(13)-VH.b(14:16)-L(15)-VL.b(16:14)|VL.a(17:6)-CL(18:7){1}-L(19)-VH.b(20:22)-L(21)-VL.b(22:20)",
    "scFV-H-scFV": "VL.a(1:3)-L(2)-VH.a(3:1)-H(4:11){2}-VH.b(5:7)-L(6)-VL.b(7:5)|VL.c(8:10)-L(9)-VH.c(10:8)-H(11:4){2}-VH.d(12:14)-L(13)-VL.d(14:12)",
    "F(ab)2":"VH.a(1:4) - CH1(2:5){1} -H(3:8){2} | VL.a(4:1) - CL(5:2){1} | VH.b(6:9) - CH1(7:10){1} - H(8:3){2} | VL.b(9:6) - CL(10:7){1}",
    "scFV":"VH.a(1:3)-L(2)-VL.a(3:1)",
    "scFV4":"VH.a(1:3)-L(2)-VL.a(3:1)-CL(4:8){1}|VL.a(5:7)-L(6)-VH.a(7:5)-CH1(8:4){1}-H(9:20){2}-CH2(10:21)-CH3(11:22)|VH.b(12:14)-L(13)-VL.b(14:12)-CL(15:19){1}|VL.b(16:18)-L(17)-VH.b(18:16)-CH1(19:15){1}-H(20:9){2}-CH2(21:10)-CH3(22:11)",
    "sdFV4":"VHH.a(1) -CH1(2:7){1} -H(3:10){2}-CH2(4:11) -CH3*@(5:12) | VHH.a(6)  -CL(7:2){1} | VHH.b(8) -CH1(9:14){1}-H(10:3){2}-CH2(11:4) -CH3*>(12:5) | VHH.b(13) -CL(14:9){1}",
    "Nanobody":"VHH.a(1)",
    "BiTE":"VH.a(1:7)-L(2)-VL.a(3:5)-L(4)-VH.b(5:3)-L(6)-VL.b(7:1)",
    "HSAbody":"VL.a(1:3)-L(2)-VH.a(3:1)-L(4)-X(5)[TYPE:FUSION, NOTE: human serum albumin]-L(6)-VH.b(7:9)-L(8)-VL.b(9:7)",
    "Cov-X-body":"X(1)[TYPE:OTHER, NOTE: pharmacophore peptide heterodimer]-VH.a(2:7)- CH1(3:8){1}-H(4:12){2}-CH2(5:11)-CH3(6:12) | VL.a(7:2)-CL(8:3){1} | X(9)[TYPE:OTHER, NOTE: pharmacophore peptide heterodimer]-VH.b(10:15)-CH1(11:16){1}-H(12:4){2}- CH2(13:5)-CH3 (14:6) | VL.b(15:10)-CL(16:11){1}",
    "Diabody":"VH.a(1:4)-L(2)-VH.b(3:6)|VL.a(4:1)-L(5)-VL.b(6:3)",
    "Miniantibody":"VH.a(1:3)-L(2)-VL.a(3:1)-H*(4:9){1}[MOD:DISULPHIDE]-X(5:10)[TYPE:ZIPPER]|VH.b(6:8)-L(7)-VL.b(8:6)-H*(9:4){1}[MOD:DISULPHIDE]-X(10:5)[TYPE:ZIPPER]",
    "scDiabody":"VH.a(1:7)-L(2)-VL.a(3:5)-L(4)-VH.b(5:3)-L(6)-VL.b(7:1)",
    "Tandem scFv-Fc":"VH.b(1:3)-L(2)-VL.b(3:1)-L(4)-VH.a(5:7)-L(6)-VL.a(7:5)-H(8:18){2}-CH2(9:19)-CH3(10:20)|VH.b(11:13)-L(12)-VL.b(13:11)-L(14)-VH.a(15:17)-L(16)-VL.a(17:15)-H(18:8){2}-CH2(19:9)-CH3(20:10)",
    "scDiabody-Fc":"VL.a(1:7)-L(2)-VL.a(3:5)-L(4)-VH.a(5:3)-L(6)-VH.a(7:1)-H(8:18){2}-CH2(9:19)-CH3(10:20)|VL.b(11:17)-L(12)-VL.b(13:15)-L(14)-VH.b(15:13)-L(16)-VH.b(17:11)-H(18:8){2}-CH2(19:9)-CH3(20:10)",
    "scDiabody-CH3":"VL.a(1:7)-L(2)-VL.a(3:5)-L(4)-VH.a(5:3)-L(6)-VH.a(7:1)-H(8:17){2}-CH3(9:18) |VL.b(10:16)-L(11)-VL.b(12:14)-L(13)-VH.b(14:12)-L(15)-VH.b(16:10)-H(17:8){2}-CH3(18:9)",
    "Diabody-CH3":"VH.b(1:11)-L(2)-VL.a(3:13){1}-H(4:9){2}-CH3(5:10)|VH.b(6:14)-L(7)-VL.a(8:16){1}-H(9:4){2}-CH3(10:5)|VL.b(11:1)-L(12)-VH.a(13:3){1}|VL.b(14:6)-L(15)-VH.a(16:8){1}",
    "DART":"VL.a(1:7)-L(2)-VH.b(3:5)-H*(4:8){3}[MOD:DISULPHIDE]|VL.b(5:3)-L(6)-VH.a(7:1)-H*(8:4){3}[MOD:DISULPHIDE]",
    "Tandem A and B": "VH.a(1:8)-L(2)-VL.b(3:10)-L(4)-VH.b(5:12)-L(6)-VL.a(7:14)|VL.a(8:1)-L(9)-VH.b(10:3)-L(11)-VL.b(12:5)-L(13)-VH.a(14:7)",
    "Intrabody":"VL.a(1:3)-L(2)-VH.a(3:1)-H(4:14){2}-CH1(5:15)-CH2(6:16)-L(7)-VH.b(8:10)-L(9)-VL.b(10:8)|VL.a(11:13)-L(12)-VH.a(13:11)-H(14:4){2}-CH1(15:5)-CH2(16:6)-L(17)-VH.b(18:20)-L(19)-VL.b(20:18)",
    "Fv-Fc":"VH.a(1:5){1}-H(2:7){2}-CH1(3:8)-CH2(4:9)| VL.a(5:1){1}|VH.b(6:10){1}-H(7:2){2}-CH1(8:3)-CH2(9:4)|VL.b(10:6){1}",
    "Triplebody":"VH.a(1:7)-CH1(2:8){1}-L(3)-VL.b(4:6)-L(5)-VH.b(6:4)|VL.a(7:1)-CL(8:2){1}-L(9)-VL.c(10:12)-L(11)-VH.c(12:10)",
    "scTriplebody":"VH.a(1:8)-CH1(2:9){2}-L(3)-VL.b(4:6)-L(5)-VH.b(6:4)-L(7)-VL.a(8:1)-CL(9:2){2}-L(10)-VL.c(11:13)-L(12)-VH.c(13:11)",
    "TriBiMinibody":"VH.a(1:3)-L(2)-VL.a(3:1)-H(4:13){2}-CH3*@(5:14)-L(6)-VH.b(7:9)-L(8)-VL.b(9:7)|VH.c(10:12)-L(11)-VL.c(12:10)-H(13:4){2}-CH3*>(14:5)",
    "LUZ-Y":"VL.a(1:4)-CL(2:5){1}-L(3)-VH.a(4:1)-CH1(5:2){1}-H(6:15){2}-CH2(7:16)-CH3(8:17)-X(9:18)[TYPE:ZIPPER]|VL.b(10:13)-CL(11:14){1}-L(12)-VH.b(13:10)-CH1(14:11){1}-H(15:6){2}-CH2(16:7)-CH3(17:8)-X(18:9)[TYPE:ZIPPER]",
    "KIH IgG-scFab":"VH.a(1:14)-CH1(2:15){1}-H(3:11){2}-CH2(4:12)-CH3*>(5:13){1}[MOD:DISULFIDE]-L(6)-VH.b(7:18)-CL*(8:19){1}[MOD:DISULFIDE] |VH.a(9:16)-CH1(10:17)-H(11:3){2}-CH2(12:4)-CH3*@(13:5){1}[MOD:DISULFIDE] |VL.a(14:1)-CL*(15:2){1}[MOD:DISULFIDE] |VL.a(16:9)-CL*(17:10)[MOD:DISULFIDE] |VL.b(18:7)-CH1(19:8){1}",
    "Dock and Lock":"VH.a(1:5)-CH1(2:6){1}-L(3)-X(4)[TYPE:FUSION, NOTE:DDD2/AD2 heterodimer]|VL.a(5:1)-CL(6:2){1}|VH.b(7:10)-CH1(8:11){1}-L(9)-X(4)[TYPE:FUSION, NOTE:DDD2/AD2 heterodimer]|VL.b(10:7)-CL(11:8){1}|VH.c(12:15)-CH1(13:16){1}-L(14)-X(4)[TYPE:FUSION, NOTE:DDD2/AD2 heterodimer]|VL.c(15:12)-CL(16:13){1}|VH.d(17:20)-CH1(18:21){1}-L(19)-X(4)[TYPE:FUSION, NOTE:DDD2/AD2 heterodimer]|VL.d(20:17)-CL(21:18){1}",
    "scFV-IgG-scFV-scFV": "VL.a(1:9)-CL(2:10){1}|VL.a(3:26)-CL(4:27){1}|VL.b(5:7)-L(6)-VH.b(7:5)-L(8)-VH.a(9:1)-CH1(10:2){1}-H(11:28){2}-CH2(12:29)-CH3(13:30)-L(14)-VH.b(15:17)-L(16)-VL.b(17:15)-L(18)-VH.c(19:21)-L(20)-VL.c(21:19)|VL.b(22:24)-L(23)-VH.b(24:22)-L(25)-VH.a(26:3)-CH1(27:4){1}-H(28:11){2}-CH2(29:12)-CH3(30:13)-L(31)-VH.b(32:34)-L(33)-VL.b(34:32)-L(35)-VH.c(36:38)-L(37)-VL.c(38:36)",
    "scFV-scFV-Fc":"VH.a(1:3)-L(2)-VL.a(3:1)-L(4)-VH.b(5:7)-L(6)-VL.b(7:5)-CH2(8:11)-CH3(9:12)-L(10)-CH2(11:8)-CH3(12:9)",
    "Trimeric Fusion Protein":"VH.a(1:6)-CH1(2:7){1}-H(3:11){2}-CH2(4:12)-CH3(5:13)|VL.a(6:1)-CL(7:2){1}|X(8:9,14)[NOTE:FUSION]-X(9:8,14)[NOTE:FUSION]-CH1(10:15){1}-H(11:3){2}-CH2(12:4)-CH3(13:5)|X(14:8,9)[NOTE:FUSION]-CL(15:10){1}",
    'scFV-X-Fc-Body':"VL.a(1:3)-L(2)-VH.a(3:1)-X(4:10){1}[TYPE:FUSION]-CH2(5:11)-CH3(6:12)|VL.b(7:9)-L(8)-VH.b(9:7)-X(10:4){1}[TYPE:FUSION]-CH2(11:5)-CH3(12:6)",
    "IgG-IgG":"VH.a(1:13)-CH1(2:14){1}-H(3:8){2}-CH2(4:9)-CH3(5:10)|VH.a(6:15)-CH1(7:16){1}-H(8:3){2}-CH2(9:4)-CH3(10:5)-L(11)-C(12)[TYPE:OPDM]|VL.a(13:1)-CL(14:2){1}|VL.a(15:6)-CL(16:7){1}|VH.b(17:28)-CH1(18:29){1}-H(19:25){2}-CH2(20:26)-CH3(21:27)-L(22)-C(12)[TYPE:OPDM]|VH.b(23:30)-CH1(24:31){1}-H(25:19){2}-CH2(26:20)-CH3(27:21)|VL.b(28:17)-CL(29:18){1}|VL.b(30:23)-CL(31:24){1}"
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
        global root
        top = tk.Toplevel(root)
        top.title = "Settings"
        top.geometry('500x400')
        tabControl = ttk.Notebook(top)

        tab1 = ttk.Frame(tabControl)
        tab2 = ttk.Frame(tabControl)

        tabControl.add(tab1, text ='Pairing Sensitivity')
        tabControl.add(tab2, text ='Colour Changer')
        tabControl.pack(expand = 1, fill ="both")
    ##Frame 1

        def Update_settings():
            global Pairing_sensitivity
            Pairing_sensitivity = sensitivity_scalebar.get()
            global Bond_thickness
            Bond_thickness = bond_thickness_scalebar.get()
            global H_Labels
            if H_label_scalebar.get() == 1:
                H_Labels = True
            else:
                H_Labels = False
            global L_Labels
            if L_label_scalebar.get() == 1:
                L_Labels = True
            else:
                L_Labels = False
            global Show_Leucine_Zippers
            if Leucine_zipper_scalebar.get() == 1:
                Show_Leucine_Zippers = True
            else:
                Show_Leucine_Zippers = False
            global Bond_Arrows
            if Bond_arrows_scalebar.get() == 1:
                Bond_Arrows = (8,10,3)
            else:
                Bond_Arrows = (0,0,0)


        settings_frame = tk.Frame(tab1, bg = "#D3D3D3")
        settings_frame.place(relx=0.05, rely = 0.05,relheight = 0.9, relwidth = 0.9)

        ttk.Label(settings_frame,text ="Pairing Sensitivity (pixels)").place(relx=0.1, rely = 0.02)
        sensitivity_scalebar = tk.Scale(settings_frame, from_=0, to=100, orient="horizontal")
        sensitivity_scalebar.set(30)
        sensitivity_scalebar.place(relx=0.1, rely = 0.10,relheight = 0.2, relwidth = 0.8)

        ttk.Label(settings_frame,text ="Bond Thickness (pixels)").place(relx=0.1, rely = 0.30)
        bond_thickness_scalebar = tk.Scale(settings_frame, from_=0, to=5, orient="horizontal")
        bond_thickness_scalebar.set(2)
        bond_thickness_scalebar.place(relx=0.1, rely = 0.37,relheight = 0.2, relwidth = 0.8)

        ttk.Label(settings_frame,text ="Bond Arrows").place(relx=0.1, rely = 0.6)
        Bond_arrows_scalebar = tk.Scale(settings_frame, from_=0, to=1, orient="horizontal")
        Bond_arrows_scalebar.set(0)
        Bond_arrows_scalebar.place(relx=0.1, rely = 0.67,relheight = 0.2, relwidth = 0.2)

        ttk.Label(settings_frame,text ="H Labels").place(relx=0.3, rely = 0.6)
        H_label_scalebar = tk.Scale(settings_frame, from_=0, to=1, orient="horizontal")
        H_label_scalebar.set(0)
        H_label_scalebar.place(relx=0.3, rely = 0.67,relheight = 0.2, relwidth = 0.2)

        ttk.Label(settings_frame,text ="L Labels").place(relx=0.5, rely = 0.6)
        L_label_scalebar = tk.Scale(settings_frame, from_=0, to=1, orient="horizontal")
        L_label_scalebar.set(0)
        L_label_scalebar.place(relx=0.5, rely = 0.67,relheight = 0.2, relwidth = 0.2)

        ttk.Label(settings_frame,text ="Leucine Zippers").place(relx=0.7, rely = 0.6)
        Leucine_zipper_scalebar = tk.Scale(settings_frame, from_=0, to=1, orient="horizontal")
        Leucine_zipper_scalebar.set(0)
        Leucine_zipper_scalebar.place(relx=0.7, rely = 0.67,relheight = 0.2, relwidth = 0.2)

        Update_settings_button= tk.Button(settings_frame, font=20, text = "Update", command =lambda: Update_settings())
        Update_settings_button.place(relx=0.333, rely = 0.90,relheight = 0.07, relwidth = 0.333)


        #ttk.Label(settings_frame,text="Hinge Labels").place(relx=0.1,rely = 0.7)

    ##Frame 2
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
        ColoursettingsLibrary.place(relx=0.1, rely = 0.1,relheight = 0.8, relwidth = 0.3)

        changerframe = tk.Frame(tab2, bg = "#D3D3D3")
        changerframe.place(relx=0.4, rely = 0.1,relheight = 0.8, relwidth = 0.49)

        colour_checker = tk.Frame(changerframe, pady = 5, bd = 5,highlightbackground="black", highlightthickness=1)
        colour_checker.place(relx=0.2, rely = 0.2,relheight = 0.2, relwidth = 0.6)

        changecolourbutton = tk.Button(changerframe, font=20, text = "Change colour", command =lambda: browse_colour())
        changecolourbutton.place(relx=0.25, rely = 0.5,relheight = 0.1, relwidth = 0.5)

        revertcolorbutton = tk.Button(changerframe, font=20, text = "Revert colour", command =lambda: revertcolor())
        revertcolorbutton.place(relx=0.25, rely = 0.7,relheight = 0.1, relwidth = 0.5)

        revertallcoloursbutton= tk.Button(changerframe, font=20, text = "Revert all colours", command =lambda: revertallcolours())
        revertallcoloursbutton.place(relx=0.25, rely = 0.9,relheight = 0.1, relwidth = 0.5)

    #Coloursettings.place(relx=0.21, rely = 0.425,relheight = 0.2, relwidth = 0.19)
    def open_manual():
        global root
        top2 = tk.Toplevel(root)
        top2.title = "README.txt"
        top2.geometry('500x400')
        manual_frame = tk.Frame(top2, bg = "#D3D3D3")
        manual_frame.place(relx=0, rely = 0,relheight = 1, relwidth = 1)
        Manual =[
        "abYdraw",
        " ",
        " ",
        "This is a programme designed to use our group's Antibody Markup Language (AbML) for describing bispecific antibody (BsAb) formats by either inputting an AbML descriptor string of a BsAb or by drawing a BsAb and outputting the its descriptor string. It is written in Python 3 and using standard packages TKinter to build the graphical interface in order to make it as accessible as possible.",
        " ",
        "Contents",
        "1. Installing and Executing",
        "2. Interface",
        "3. AbML",
        "4. Inputting AbML",
        "5. Obtaining AbML",
        "6. Formats Library",
        "7. Saving and exporting",
        "8. Settings",
        " ",
        "1. Installing and Executing",
        "abYdraw may be downloaded and run as an executable Python script. Therefore it requires at least Python 3.8 to run and can be executed as such:",
        "",
        "python3 abYdraw.py"
        " ",
        "2. Interface",
        "The programme interface includes six points of reference, four of which in a column on the left hand side and two more on the right hand side.",
        "Starting with the left hand column, the first is the Domain palette which has buttons necessary for drawing antibody domains, secondly a library of commonly used bispecific antibody AbML expressions, thirdly the input box for AbML expressions and a buttonpad that will render antibody schematics or output AbML to the textbox.",
        "On the right hand side, the most prominent feature is the canvas for drawing and rendering antibody schematics and underneath there are two buttons which are involved in exporting the schematic.",
        " ",
        "3. AbML",
        "Our language was derived from existing macromolecule descriptor languages but we have compensated for their limitations and made AbML simple whilst conveying as much useful information as possible. Strings are split into chains, which are then split into domains. Each domain type has its own symbol and each domain unit also carries additional information including: modification types; the specificity of the variable region (if applicable); a number label assigned to the domain and the",
        "number label assigned to the domain it interacts with; the number of disulphide bonds between the two interacting domains and comments outlining additional information not covered by the language of types: TYPE; NOTE; MOD; ANTI and LENGTH. Full descriptions of AbML can be found on the language guide sheet included in the Repository.",
        " ",
        "4. Inputting AbML",
        "AbML descriptor strings my be inputted in the entry box or opened in the 'File>Open' menu and then by clicking 'Get Structure', a schematic of that antibody will render in the canvas. Schematics are drawn in colour-coded fashion depending on any specificities given in the descriptor chain. Domains are connected by different kinds of linkers which are also colour-coded depending on their type. Any comments given in the descriptor string are also displayed beside the schematic. Labels on",
        "the schematic may be toggled on and off using the Labels key on the Domain Palette.",
        " ",
        "5. Obtaining AbML",
        "To draw an antibody, you must insert domains onto the canvas and arrange them so the programme recognises it as an antibody. Your tools for adding domains to the canvas are in the Domain Palette which contains all of the domain types, modifications, specificities and comment types as described in the AbML guidesheet as well as some options. Selecting a button on the palette will cause it to flash red to indicate it has been switched on. Only one domain or connector type may be switched",
        "on at a time, but you may choose any combination of modifications and specificities to accompany your selection. Specificity types are only applicable to VH or VL domains but when drawing other domains, a default 'a' value is set. Once selected you will notice the cursor will change from arrow to '+' sign. This means you can left-click to insert your chosen domain type onto the canvas at the location you have clicked. If no domain types are selected but a modification or specificity",
        "are, then by clicking on a domain on the canvas, you may replace its current specificity and modification to those you have selected.",
        "Domains may be connected with the connector options in the first column of the palette. When a connector type is selected it will become highlighted in red and the user must click and drag the bond from its starting position to its end, making sure each end is inside the boundaries two domains it links. Bonds are unidirectional and start from N-terminus to C-terminus.",
        "Comments may be added by selecting a comment type, which will highlight the button just pressed and the 'Comment' button and then inputting the comment into the entry box beneath the palette. Comments may then be drawn on the domains they are applicable to. To disable further commenting, ensure the comment type and 'Comment' buttons are no longer highlighted.",
        "If no domain types or modifications are selected then domains, linkers and comments may be rearranged by clicking and dragging them. You must ensure that interacting domains are positioned close together at roughly the same level and that VH and VL domains face each other to complete their antigen binding domains. Their orientations can be changed by right-clicking the relevant domains. Furthermore features may be deleted by selecting the 'Delete' button on the palette and",
        "then selecting what to delete. To remove all features, click 'Clear All' and the canvas will be made blank.",
        "Once domains and connectors are arranged, click the 'Get AbML' button to generate the AbML descriptor string for this sequence. Once this has been generated it will appear in the input box. You may then click 'Get structure' again to re-render the schematic with abYdraw. Alternatively, the 'Tidy' button performs both steps of this operation. Once rendered, an image may be altered by adding or removing domains. By clicking 'Get Sequence' or 'Tidy', you will obtain a new expression",
        "for rendering.",
        " ",
        "6. AbML Formats Library",
        "To assist users, the programme has a library of BsAb formats available which can be scrolled through and selected. This will give the schematic and AbML expression for this format that can be used as a starting point to make new expressions and schematics that are relevant to the user.",
        " ",
        "7. Exporting and Saving",
        "Exporting the canvas image as .eps file can be done by 'Export EPS' and using the file directory to save the image. Alternatively the AbML may be saved as a text file by clicking the 'File>Save' option in the menu.",
        "Finally you may export a Template File from the expression in the entry box which is a format of notating important BsAb residues. The programme cannot locate these residues but it can identify the features that are in the BsAb that users may want to include in the Template Files.",
        " ",
        "8. Settings",
        "Users may change aspects about the rendering and pairing of their schematic. By opening the 'File>Settings' window, a menu with two tabs will appear. The first tab has settings regarding pairing sensitivity of drawn schematics, bond thickness, directional arrows, Hinge and Linker labels which can be set by using the appropriate sliders. For pairing sensitivity the scale is 0-80 pixels and bond width are between 1-5 pixels. Other binary settings are on sliders 0-1.",
        "To update the settings,users must press the update button and re-render their schematic to see their new schematic. The second tab is the colour-coding menu which with a list of domain types. When a domain type is selected the current colour of assigned to that domain will appear. 'Change colour' allows users to assign a new colours to that domain type using the colour palette of the operating system, but these may be reverted by 'Revert colour' or 'Revert all colours'",
        " ",
        ]
        manual_textbox = tk.Text(manual_frame,wrap="word", font=20)
        manual_textbox.place(relx=0, rely = 0, relheight = 1, relwidth = 1)
        for i in range(len(Manual)):
            manual_textbox.insert("end",str(Manual[i]+"\n"))
    ###Results canvas

    lower_frame = tk.Frame(root, bg = '#80c1ff', bd=5)
    lower_frame.place(relx=0.45, rely=0.015, relwidth=0.55,relheight=0.93)
    lower_frame2 = tk.Frame(lower_frame, width =  700+200, height = 700+300, bg = '#80c1ff', bd=5)
    #lower_frame2.place(relwidth=1,relheight=1)
    lower_frame2.pack(expand=True,fill="both")
    lower_canvas = tk.Canvas(lower_frame2,width=700+200,height=700+300, scrollregion=(0,0,700+200,700+300))
    yscrollbar = tk.Scrollbar(lower_frame2, orient="vertical")
    yscrollbar.config(command=lower_canvas.yview)
    yscrollbar.pack(side="right",fill="y")
    xscrollbar = tk.Scrollbar(lower_frame2, orient="horizontal")
    xscrollbar.config(command=lower_canvas.xview)
    xscrollbar.pack(side="bottom",fill="x")
    lower_canvas.config(yscrollcommand=yscrollbar.set)
    lower_canvas.config(xscrollcommand=xscrollbar.set)
    #lower_canvas.place(relheight=1,relwidth=1)
    lower_canvas.pack(expand=True,fill="both")


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
    lower_canvas.config(cursor = "fleur")
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.bind("<Button-2>", mm.change_orientation)
    lower_canvas.bind("<Button-3>", mm.change_orientation)
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
    template_file_button = tk.Button(export_frame, text = "Export template file", bg = "#CFCFCF", font=20, command=lambda: Get_Template_File(lower_canvas))
    template_file_button.place(relx=0, rely=0, relwidth=0.5,relheight=1)
    Image_file_button = tk.Button(export_frame, text = "Export EPS", bg = "#CFCFCF", font=20, command=lambda: save_as_png(lower_canvas))
    Image_file_button.place(relx=0.5, rely=0, relwidth=0.5,relheight=1)
    #img = str("AbYdraw_icon.png")
    #root.iconphoto(False, tk.PhotoImage(file=img))



    tite = root.title('abYdraw')
    menubar = tk.Menu(root)
    filemenu = tk.Menu(menubar, tearoff=0)
    filemenu.add_command(label="Settings", command=open_settings)
    filemenu.add_command(label="Open", command=lambda: browseFiles())
    filemenu.add_command(label="Save", command=lambda: save_txt_file())
    filemenu.add_command(label="Export template file", command=lambda: Get_Template_File(lower_canvas))
    filemenu.add_command(label="Export EPS", command=lambda: save_as_png(lower_canvas))

    filemenu.add_separator()

    filemenu.add_command(label="Exit", command=root.quit)
    menubar.add_cascade(label="File", menu=filemenu)

    editmenu = tk.Menu(menubar, tearoff=0)
    editmenu.add_command(label="Undo", command=lambda: undo())
    editmenu.add_command(label="Redo", command=lambda: redo())
    menubar.add_cascade(label="Edit", menu=editmenu)

    helpmenu = tk.Menu(menubar, tearoff=0)
    helpmenu.add_command(label="Manual", command=open_manual)
    menubar.add_cascade(label="Help", menu=helpmenu)

    root.config(menu=menubar)


    ##Choose between GUI and CLI

    root.mainloop()
else:
    # Create the parser
    my_parser = argparse.ArgumentParser(prog="abYdraw",
                                        usage='python %(prog)s [options]',
                                        description='Draws antibody formats')

    my_parser.add_argument('-f','--file', type=str,help='the path to plaintext file with AbML expression')
    my_parser.add_argument('-i','--input', type=str,help='string of AbML input')
    my_parser.add_argument('-o','--output', type=str,help='string of image output name (default "abYdraw_export")')
    my_parser.add_argument('-s','--show',nargs='?', type=int,help='Show image window 0-1 (default 0)')
    my_parser.add_argument('-e','--format',type=str ,help='specify image format eps, png, jepg (default eps)')
    my_parser.add_argument('-p','--image',nargs='?', type=int,help='Save image file 0-1 (default 1)')
    my_parser.add_argument('-t','--template',nargs='?', type=int,help='Save template file 0-1 (default 0)')
    my_parser.add_argument('-l','--labels',nargs='?', type=int,help='Toggle domain labels 0-1 (default 1)')
    my_parser.add_argument('-j','--hinge', nargs='?',type=int,help='Toggle hinge labels 0-1 (default 0)')
    my_parser.add_argument('-k','--linker', nargs='?', type=int,help='Toggle linker labels 0-1 (default 0)')
    my_parser.add_argument('-a','--arrows', nargs='?',type=int,help='Toggle bond direction arrows 0-1 (default 0)')
    my_parser.add_argument('-b','--thickness', type=int,help='Set bond thickness 1-5 (default 2)')

    my_parser.print_help()
    args = my_parser.parse_args()
    print(args)
    input_path = args.file
    input_string = args.input
    if input_string is not None and input_path is None:
        entry = input_string
    if input_path is not None and input_string is None:
        entry = Get_input(input_path)
    elif input_string is not None and input_path is not None:
        entry = input_string
    elif input_string is None and input_path is None:
        print('No input was given. Exiting programme')
        sys.exit()
    output_name = args.output
    if output_name is None:
        output_name = "AbYdraw_export"

    Label_lock = args.labels
    if Label_lock is None or Label_lock != 0:
        Label_lock = 1
    elif Label_lock == 0:
        Label_lock = 0

    H_Labels = args.hinge
    if H_Labels is None or H_Labels != 0:
        H_Labels = 0
    elif H_Labels == 0:
        H_Labels = 1

    L_Labels = args.linker
    if L_Labels is None or L_Labels != 0:
        L_Labels = 0
    elif L_Labels == 0:
        L_Labels = 1


    Bond_Arrows = args.arrows
    if Bond_Arrows is None:
        Bond_Arrows = (0,0,0)
    elif Bond_Arrows is not None:
        Bond_Arrows = (8,10,3)

    Bond_thickness = args.thickness
    if Bond_thickness == None:
        Bond_thickness = 2

    show = args.show
    if show is None:
        show =0
    elif show is not None:
        show = 1

    template = args.template
    if template is None:
        template = 0
    elif template is not None:
        template = 1

    save = args.image
    if save is None or save != 0:
        save = 1
    elif save is not None:
        save = 0
    format = args.format
    if format is None:
        format = "eps"

    CLI = True
    all_buttons = []
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
    canvas_polygons = {}
    canvas_labels   = {}
    temp_label      = {}
    TYPE_labels     = {}
    NOTE_labels     = {}
    ANTI_labels     = {}
    MOD_labels      = {}
    LENGTH_labels   = {}

    specificity_colours = [a_heavy_colour, a_light_colour,b_heavy_colour, b_light_colour,c_heavy_colour, c_light_colour,d_heavy_colour, d_light_colour,e_heavy_colour, e_light_colour,f_heavy_colour, f_light_colour,g_heavy_colour, g_light_colour,h_heavy_colour, h_light_colour,generic_heavy_colour,generic_light_colour,X_colour,C_colour]
    bond_colours = [disulphide_colour,bond_colour,hinge_colour,linker_colour]


    def CLI_save_png(format, output,canvas):
        fileName = str(output)
        eps = canvas.postscript(file=(fileName+"."+format), colormode='color', height= 1000, width = 900)
        print(str(fileName+"."+format))
        abs_path = str(fileName+"."+format)
        if "eps" not in format:
            img = Image.open(str(abs_path))
            img.load(scale=10)
            if img.mode in ('P', '1'):
                img = img.convert("RGB")
            TARGET_BOUNDS = (1024, 1024)
            ratio = min(TARGET_BOUNDS[0] / img.size[0], TARGET_BOUNDS[1] / img.size[1])
            new_size = (int(img.size[0] * ratio), int(img.size[1] * ratio))
            img = img.resize(new_size, Image.ANTIALIAS)
            if "png" in format:
                img.save(abs_path, "png")
            elif "jpeg" in format:
                img.save(abs_path, "png")

    def CLI_function(input_string,output_name,Label_lock, H_Labels, L_Labels, Bond_Arrows, Bond_thickness, show, template, save, format):
        ###Window###

        HEIGHT =  900
        WIDTH  = 1000
        root.geometry=("900x1000")
        lower_frame = tk.Frame(root, width =  700+200, height = 700+300, bg = '#80c1ff', bd=5)
        lower_frame.pack(expand=True,fill="both")
        lower_canvas = tk.Canvas(lower_frame,width=700+200,height=700+300, scrollregion=(0,0,700+200,700+300))
        yscrollbar = tk.Scrollbar(lower_frame, orient="vertical")
        yscrollbar.config(command=lower_canvas.yview)
        yscrollbar.pack(side="right",fill="y")
        xscrollbar = tk.Scrollbar(lower_frame, orient="horizontal")
        xscrollbar.config(command=lower_canvas.xview)
        xscrollbar.pack(side="bottom",fill="x")
        lower_canvas.config(yscrollcommand=yscrollbar.set)
        lower_canvas.config(xscrollcommand=xscrollbar.set)
        lower_canvas.pack(expand=True,fill="both")
        Bond_thickness = Bond_thickness
        text_to_image = True
        split_chains = Get_dictionaries(entry)
        coordinates  = Check_interactions(split_chains, lower_canvas)
        rendering = render(coordinates, lower_canvas,True)
        if format != "png" and format != "jpeg" and format !="eps":
            message = "Please specify a correct format"
            raise_error(message,lower_canvas)
        if save == 1:
            CLI_save_png(format, output_name,lower_canvas)
        if template == 1:
            Get_Template_File(lower_canvas)
        if show == 1:
            root.mainloop()

    CLI_function(input_string,output_name,Label_lock, H_Labels, L_Labels, Bond_Arrows, Bond_thickness, show, template, save, format)
