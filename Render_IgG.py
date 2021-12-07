#!/usr/bin/python
import re
import sys
import os
from PIL import Image
import tkinter as tk

HEIGHT = 1200
WIDTH  = 1600

####If numbers cut out, need to be doubled
######################################
def Get_input(x):
    input = ""
    with open(x,"r") as f:
        for line in f:
            input += line

    return input
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
    if "|" in y:
        splitx  = y.split("|")
    else:
        splitx = [y]
    if len(splitx) == 4:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a']
    elif len(splitx) == 3:
        chains  = ['VH.b', 'VL.b','VH.a']
    elif len(splitx) == 2:
        chains  = ['VH.b','VH.a']
    elif len(splitx) == 1:
        chains  = ['VH.a']
    elif len(splitx) ==5:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1']
    elif len(splitx) ==6:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1', 'fragment2']
    elif len(splitx) ==7:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1', 'fragment2', 'fragment3']
    elif len(splitx) ==8:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1', 'fragment2', 'fragment3', 'fragment4']

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
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                elif chain[1] != "L":
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                elif chain[1] == "L":
                    if "VH" in chain[2]:
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                    elif "VL" in chain[2]:
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))

            elif "VL" in chain[0]:
                try:
                    if  chain[1] != "L":
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                    elif chain[1] == "L":
                        if len(chain) > 2:
                            if "VH" in chain[2]:
                                dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                            elif "VL" in chain[2]:
                                dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                except IndexError:
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))

            elif "X" in chain[0]:
                 dict = str(re.sub("\.|\+|\_","",str(chains[i])))
        elif i >= 4:
            dict_number = i-3
            dict = str("fragment"+str(dict_number))


        for j in range(len(chain)):
            domain   =  str(chain[j])
            domain = str(re.sub("\[.*\]|\(.*\)|\{.*\}|\[|\'|\]|\.|\*","", str(domain)))

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
            locationstr = locationstr.split(",")
            if domain != "Linker":
                for i in locationstr:
                    location.append(int(i))
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

    print(VHa)
    print(VLa)
    print(VHb)
    print(VLb)
    print(fragment1)
    ###checker###
    VHa_keyslist = list(VHa.keys())
    VLa_keyslist = list(VLa.keys())
    VHb_keyslist = list(VHb.keys())
    VLb_keyslist = list(VLb.keys())
    VHa_checked = {}
    VHb_checked = {}
    VLa_checked = {}
    VLb_checked = {}
    ##check VH chains interact

    chains = [VHa_keyslist,VLa_keyslist,VHb_keyslist,VLb_keyslist]
    dicts = [VHa,VLa,VHb,VLb]

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
                                    break
                                elif current_interactor == interactor and ("H[" in str(chains[i][j]) and "H[" in str(chains[a][b])) and ("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) and VHa_VHb_found == False:
                                    VHa_VHb_found = True
                                    VHa_checked,VHb_checked= dicts[i],dicts[a]
                                    break
                                elif current_interactor == interactor and ("H[" not in str(chains[i][j]) and "H[" not in str(chains[a][b])) and ("X[" in str(chains[i][j]) and "X[" in str(chains[a][b])) and VHa_VHb_found == False :
                                    VHa_VHb_found = True
                                    VHa_checked,VHb_checked= dicts[i],dicts[a]
                                    break
                except IndexError:
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
                                    elif i %2 != 0 and VLb_match == False:
                                        VLb_checked = dicts[a]
                                        VLb_match = True

                except IndexError:
                    continue
    elif chain_count ==2:
        VHa_checked = VHa
        VHb_checked = VHb
    elif chain_count ==1:
        VHa_checked = VHa
    return(VHa_checked,VLa_checked,VHb_checked,VLb_checked,chain_count,fragment1,fragment2,fragment3,fragment4)

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
    VHa_chain_master       = VHa_chain.copy()
    VLa_chain_master       = VLa_chain.copy()
    VHb_chain_master       = VHb_chain.copy()
    VLb_chain_master       = VLb_chain.copy()
    chain_count     = chains_list[4]
    fragment1       = chains_list[5]
    fragment2       = chains_list[6]
    fragment3       = chains_list[7]
    fragment4       = chains_list[8]
    All_positions_and_chains    ={}
    extra_disulphide_bridges    ={}
    H_disulphide_coordinates    ={}
    completed_disulphidebridges=[]
    Notes           = []
    Notes_positions = []
    noted_Hbonds=False
    H_coordinatey   = []


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
        previous_linker = 0

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
                    ADCs.append(fragments[x][i])

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
                x_cleaned_up = re.sub("\[.*\]","", x)

                return(str(x_cleaned_up))


            if "a" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "d" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_a.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_a.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_a.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_a.append(cleaned_up)
            elif "b" in str(fragments[x]) and "a" not in str(fragments[x]) and "c" not in str(fragments[x]) and "d" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_b.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_b.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_b.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_b.append(cleaned_up)
            elif "c" in str(fragments[x]) and "a" not in str(fragments[x]) and "b" not in str(fragments[x]) and "d" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_c.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_c.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_c.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_c.append(cleaned_up)
            elif "d" in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "a" not in str(fragments[x]):
                for i in range(len(fragments[x])):
                    if "H" in str(fragments[x][i]):
                        coordinates_list_heavy_d.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_heavy_d.append(cleaned_up)
                    elif "L" in str(fragments[x][i]):
                        coordinates_list_light_d.append(dictionary.get(fragments[x][i]))
                        cleaned_up = fragment_cleanup(fragments[x][i])
                        names_list_light_d.append(cleaned_up)
            elif "a" not in str(fragments[x]) and "b" not in str(fragments[x]) and "c" not in str(fragments[x]) and "d" not in str(fragments[x]):
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
                    else:
                        if "H" in str(fragments[x][i]):
                            coordinates_list_heavy_a.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_heavy_a.append(cleaned_up)
                        elif "L" in str(fragments[x][i]):
                            coordinates_list_light_a.append(dictionary.get(fragments[x][i]))
                            cleaned_up = fragment_cleanup(fragments[x][i])
                            names_list_light_a.append(cleaned_up)
        print("OKILY DOKILY", coordinates_list_heavy_a)
        print("DOOBIE DOOBIE", coordinates_list_heavy_a)
        return(coordinates_list_heavy_a,coordinates_list_light_a,coordinates_list_heavy_b,coordinates_list_light_b,coordinates_list_heavy_c,coordinates_list_light_c,coordinates_list_heavy_d,coordinates_list_light_d,names_list_heavy_a,names_list_light_a,names_list_heavy_b,names_list_light_b,names_list_heavy_c,names_list_light_c,names_list_heavy_d,names_list_light_d)


    def innie_or_outie(chain,VHa_chain,VHb_chain,VLa_chain,VLb_chain,Build_in,Build_out,fragment1,fragment2,fragment3,fragment4, righthanded):
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
                    print("CUM ERE LUV 2")
                    check_VHa_VHb_interactor = False
                    second_comp_keyslist = list(second_comp.keys())
                    for i in range(len(second_comp_keyslist)):
                        current_interactor = second_comp.get(second_comp_keyslist[i])[0][0]
                        if interactor == current_interactor:
                            check_VHa_VHb_interactor = True
                    if check_VHa_VHb_interactor == True:
                        print("WHY OH WHY OH WHY")
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
                                    print("OH MY GOD WE FOUND IT")
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
                                    print("WHat a fucktard")
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
                                                    print("UMM WTF")
                                                    innie_or_outie_list.append("innie")
                                            else:
                                                print("OOOH FFS")
                                                innie_or_outie_list.append("outie")
                                        except IndexError:
                                            print("oopsie")
                                            innie_or_outie_list.append("outie")
                                    elif previous_domain_interaction == True:
                                        if innie_or_outie_list[-2] == "outie":
                                            print("OUT THOT")
                                            innie_or_outie_list.append("innie")
                                        elif innie_or_outie_list[-2] == "innie":
                                            print("IN THOT")
                                            innie_or_outie_list.append("outie")
                                    elif innie_or_outie_list[-2] == "outie":
                                        print("OUT THOT")
                                        innie_or_outie_list.append("outie")
                                    elif innie_or_outie_list[-2] == "innie":
                                        print("IN THOT")
                                        innie_or_outie_list.append("innie")

        def in_or_out_light(n,chain,first_comp,Build_in,Build_out,Light_chain_check):
            check_VHa_VLa_interactor = False
            chain_keyslist = list(chain.keys())
            interactor = chain.get(chain_keyslist[n])[0][1]
            first_comp_keyslist = list(first_comp.keys())
            for i in range(len(first_comp_keyslist)):
                current_interactor = first_comp.get(first_comp_keyslist[i])[0]
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
                else:
                    Light_chain_check = True
                    default = "innie"
                if n == 0:
                    try:
                        #print("OKILY DOKIE")
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
                                        print("LATERS LUV")
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
                    else:
                        innie_or_outie_list.append("innie")


            elif "V" not in keyslist[n] and "H[" not in keyslist[n]:
                "YEAH RIGHT"
                innie_or_outie_list.append("constant")
            elif "V" not in keyslist[n] and "H[" in keyslist[n]:
                "FFS"
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
        innie_or_outie_list=[]
        interaction_counter = 0
        keyslist = list(dictionary.keys())
        H_count = 0


        if chain_count >= 4 :
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
            if "H[" in str(dictionary) or "X[" in str(dictionary):
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
                if tangle_found == False :
                    if "X" in str(dictionary):
                        if "LEUCINE" in str(dictionary):
                            slant = True
                        else:
                            slant = False
                    elif "H[" in str(dictionary):
                        slant = True
                else:
                    slant = False
            else:
                slant = False


            keyslista = list(VHa_chain.keys())
            keyslistb = list(VHb_chain.keys())
            try:
                if dictionary == VHb_chain and VHa_chain.get(keyslista[0])[0][0] == (VHb_chain.get(keyslistb[0])[0][1]) and VHa_chain.get(keyslista[0])[0][1] == (VHb_chain.get(keyslistb[0])[0][0]):
                    print("Dear god 1")
                    Build_out = True
                    Build_in = False
                elif dictionary == VHa_chain and VHa_chain.get(keyslista[0])[0][0] == (VHb_chain.get(keyslistb[0])[1]) and VHa_chain.get(keyslista[0])[0][1] == (VHb_chain.get(keyslistb[0])[0]):
                    print("Dear god 2")
                    Build_out = True
                    Build_in = False
                else:
                    print(dictionary)
                    try:
                        if dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[2])[0][1]) and "H[" in keyslist[3]:
                            Build_in = True
                            Build_out = False
                        elif dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[2])[0][1]) and dictionary.get(keyslist[4])[0][0] == (dictionary.get(keyslist[6])[0][1]):
                            print("Dear god 3")
                            Build_out = True
                            Build_in = False
                        elif dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[6])[0][1]) and dictionary.get(keyslist[2])[0][0] == (dictionary.get(keyslist[4])[0][1]):
                            print("Dear god 4")
                            Build_in = False
                            Build_out = True
                        else:
                            print("Dear god 5")
                            Build_in = True
                            Build_out = False
                    except IndexError:
                        print("Dear god 6")
                        Build_in = True
                        Build_out = False
            except IndexError:
                print("Dear god 7")
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
        if dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2:
            Light_chain_check = True
        if dictionary == fragment1 or dictionary == fragment2:

            if starty < 200:
                slant = True
            else:
                slant = False
            if startx < 350:
                righthanded=False
            elif startx>350:
                righthanded=True
        if dictionary == VHa_chain or dictionary == VHb_chain:
            if (len(VHa_chain) == len(VHb_chain)):
                equal_chain_lengths = True
            else:
                equal_chain_lengths = False
        else:
            equal_chain_lengths = True
        innie_or_outie_list = innie_or_outie(dictionary, VHa_chain_master,VHb_chain_master,VLa_chain_master,VLb_chain_master,Build_in,Build_out, fragment1,fragment2,fragment3,fragment4,righthanded)
        print(innie_or_outie_list)

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

            Build_up=False
            Build_down=True
            Domain_name = str(re.sub("\@|\>|\<|\[.*\]","",str(keyslist[i])))


            if dictionary.get(keyslist[i])[2] != "":
                note_label = str(dictionary.get(keyslist[i])[2])
                Notes.append(Domain_name+" "+note_label)
                if len(Notes_positions) == 0:
                    Notes_positions.append([200,600])
                elif len(Notes_positions) > 0:
                    XY = (Notes_positions[-1][1])+20
                    Notes_positions.append([200,XY])

            print(mod_label, mod)

            if dictionary == VHa_chain or dictionary == VHb_chain or dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2:
            #    pass
            #elif dictionary != VHa_1_test and dictionary != VHb_1_test and dictionary != VLa_1_test and dictionary != VLb_1_test:
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
                print("checkpoint1")
                getcoordinates = domainmaker(startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H)

            elif i > 0:
                previous_domain = keyslist[i-1]
                previous_chain  = chain[i-1]
                if dictionary.get(keyslist[i])[0] != ['']:
                    checker = dictionary.get(keyslist[i])[0]
                else:
                    checker = 0
                if "H[" not in keyslist[i-1] and "Linker[" not in keyslist[i-1] and "X" not in keyslist[i-1]:
                    previous_number = (dictionary.get(previous_domain)[0])+1
                elif "X" in keyslist[i-1] and "Linker[" not in keyslist[i]:
                    print("checkpointX")
                    previous_number = (dictionary.get(keyslist[i])[0])
                    dictionary[keyslist[i-1]] = [previous_number,0]
                    previous_domain = (keyslist[i])
                elif "X" in keyslist[i-1] and chain_count > 2 and "H[" not in str(dictionary) and Light_chain_check == False:
                    slant=False



                if "H[" in keyslist[i-1]:

                    previous_domain = keyslist[i-2]
                    previous_chain  = chain[i-2]
                    print(dictionary.get(previous_domain))
                    previous_number = int(dictionary.get(previous_domain)[0])+2


                if dictionary == VLa_chain or dictionary == VLb_chain:
                    if "CL[" in keyslist[i-1]:
                        slant = False

                #print(keyslist[i], dictionary.get(keyslist[i])[0], previous_number, dictionary.get(keyslist[i])[0] , (dictionary.get(previous_domain)[1]))

                if checker in All_positions_and_chains:
                    startx = All_positions_and_chains.get(checker)[0][0]
                    starty = All_positions_and_chains.get(checker)[0][1]
                    print(startx)
                    print(startx)
                    getcoordinates = domainmaker(startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H)

                if "H[" in keyslist[i]:
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
                                getcoordinates = domainmaker((previous_chain[6])-50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant == True and righthanded == False:
                                getcoordinates = domainmaker((previous_chain[6])+50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant==False:
                                getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            Extra_bond=True
                            Build_up=True
                            Build_down=False


                        elif dictionary.get(keyslist[i])[1] != (dictionary.get(keyslist[i-1])[1]+2) or dictionary == VHb_chain:
                            print("checkpoint3")


                            if righthanded == True and slant==True:
                                getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False and slant==True:
                                getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant==False:
                                getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)



                elif "H[" in keyslist[i-1]  and "X" in keyslist[i] :
                    print("checkpoint4")
                    if righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)



                elif "X" in keyslist[i] and "H[" not in keyslist[i-1]:
                    print("checkpoint5")
                    if "Linker[" in keyslist[i-1]:
                        previous_chain = chain[i-2]
                    if "C" in keyslist[i-1]:
                        if Build_in == True:
                            Build_in = False
                            Build_out = True
                    if chain_count == 2 and  mod !="Leucine":
                        if Build_in == True:
                            Build_in = False
                            Build_out = True
                        if righthanded == False and slant == True:
                            getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and slant == True:
                            getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == False:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    elif chain_count !=2 and mod !="Leucine" and "X" not in keyslist[i-1]:
                        print("DOCCY WHOO",innie_or_outie_list[i-2])
                        if innie_or_outie_list[i-2] == "innie" and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif innie_or_outie_list[i-2] == "outie" and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif innie_or_outie_list[i-2] == "innie" and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif innie_or_outie_list[i-2] == "outie" and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                        elif righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker((previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker((previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker((previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker((previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_in = True
                        Build_out = False
                    elif chain_count !=2 and mod !="Leucine" and "X" in keyslist[i-1]:
                        if righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker((previous_chain[0]+30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker((previous_chain[0]-30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker((previous_chain[0]-30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker((previous_chain[0]+30),(previous_chain[1]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_up = True
                        Build_in = True
                        Build_out = False
                    elif mod =="Leucine":
                        if righthanded == False :
                            getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True :
                            getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == False:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)


                elif "H[" in keyslist[i-1]  and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    print("checkpoint6")
                    previous_H = True
                    slant=False
                    if chain_count <=2:
                        if dictionary.get(keyslist[2]) != [''] and "Linker[" in keyslist[2]:
                            if dictionary.get(keyslist[0])[0] != (dictionary.get(keyslist[2])[1]) :
                                getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                        else:
                            if righthanded == True :
                                getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False:
                                getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    elif H_count ==1 :
                        if righthanded == True :
                            getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif H_count==2 :
                        if righthanded == True :
                            getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    Build_in = False
                    Build_out = True

                elif "Linker[" not in keyslist[i-1] and "X" and len(dictionary.get(keyslist[i-1])) ==1:
                    print("checkpoint7")
                    if "Linker[" in keyslist[i]:
                        pass
                        print("checkpoint7.5")
                    elif slant == True and righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    else:
                        getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    Build_in  = False
                    Build_out = True

                elif "H[" not in keyslist[i] and "Linker[" not in keyslist[i-1] and "Linker[" not in keyslist[i] and "X" not in keyslist[i] and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-1])[1]+2):
                    print("checkpoint8")
                    if chain_count == 2:
                        if dictionary == VHa_chain:
                            if slant == True and righthanded == True:
                                getcoordinates = domainmaker((previous_chain[6])-50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant == True and righthanded == False:
                                getcoordinates = domainmaker((previous_chain[6])+50,(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            else:
                                getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+115),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            Build_up=True
                            Build_down=False
                        elif dictionary == VHb_chain:
                                if slant == True and righthanded == True:
                                    getcoordinates = domainmaker((previous_chain[6])-50,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif slant == True and righthanded == False:
                                    getcoordinates = domainmaker((previous_chain[6])+50,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                else:
                                    getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

                    #Build_in  = False
                    #Build_out = True


                elif "H[" not in keyslist[i] and "Linker[" not in keyslist[i-1] and "Linker[" not in keyslist[i] and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    print("checkpoint9")
                    if slant == True and righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    else:
                        getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)


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
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_in  = False
                        Build_out = True
                    elif len(dictionary.get(keyslist[i-2])) ==1 and "X" not in keyslist[i-2]:
                        print("checkpoint11")
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
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
                        if chain_count == 1 and to_joinx <= 350 and dictionary.get(keyslist[i])[1] != (dictionary.get(keyslist[i-2])[0]) and len(dictionary) >=8:
                            righthanded = True
                            Build_in = False
                            Build_out = True
                        elif chain_count == 1 and to_joinx >= 350 and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-2])[0]):
                            Build_in = False
                            Build_out = True
                        if to_join_righthanded == False and to_join_direction == 'outie':
                            print("BS1")
                            getcoordinates = domainmaker(to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == False and to_join_direction == 'innie':
                            print("BS2")
                            getcoordinates = domainmaker(to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == False and to_join_direction == 'constant':
                            print("BS2")
                            getcoordinates = domainmaker(to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == True and to_join_direction == 'outie':
                            print("BS1")
                            getcoordinates = domainmaker(to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == True and to_join_direction == 'innie':
                            print("BS2")
                            getcoordinates = domainmaker(to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif to_join_righthanded == True and to_join_direction == 'constant':
                            print("BS2")
                            getcoordinates = domainmaker(to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
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
                    elif "X" in keyslist[i-2]  :
                        print("checkpoint14")
                        if chain_count ==1:
                            if Build_in == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif righthanded == True:
                                    getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif Build_out == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif righthanded == True:
                                    getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1]+100), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif chain_count > 1:
                            if chain_count == 2:
                                slant = False
                            if righthanded == False and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker((previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker((previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker((previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker((previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            if chain_count==2:
                                Build_out = True
                                Build_in = False
                            else:
                                Build_in = True
                                Build_out = False
##Build up
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == (dictionary.get(keyslist[0])[1]):
                        print("checkpoint15")
                        getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_up=True
                        Build_down=False
                    elif chain_count == 2 and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and  dictionary.get(keyslist[i])[1]+1 == dictionary.get(keyslist[i-2])[1]:
                        print("checkpoint16")
                        if dictionary == VHa_chain:
                            if slant==True and righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0])-45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif  slant==True and righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0])+45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif slant==False:
                                getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            Build_up=True
                            Build_down=False
                        elif dictionary == VHb_chain:
                            try:
                                if "Linker[" in keyslist[i+1]:
                                    if slant==True and righthanded == True:
                                        getcoordinates = domainmaker((previous_chain[0])+45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                    elif  slant==True and righthanded == True:
                                        getcoordinates = domainmaker((previous_chain[0])-45,(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                    elif slant==False:
                                        getcoordinates = domainmaker((previous_chain[0]),(previous_chain[1])-95, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                    Build_up=True
                                    Build_down=False
                                else:
                                    getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            except IndexError:
                                getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7])+20, righthanded,slant,V,direction,X,mod,interaction,previous_H)

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
                                getcoordinates = domainmaker((previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif Build_out == True:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0]-60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]+60),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)

##Build down
                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                        print("checkpoint18")
                        if "V" in keyslist[i]:
                            in_out_counter +=1
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)

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
                                Notes_positions.append([200,600])
                            elif len(Notes_positions) > 0:
                                XY = (Notes_positions[-1][1])+20
                                Notes_positions.append([200,XY])






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
                if "H[" in keyslist[i-1] and i+1 != len(keyslist):
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
                    if righthanded==False and slant == True:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]
                        arc_topy  = top_bond[1]
                        arcbottomx= bottom_bond[0]+50
                        arcbottomy= bottom_bond[1]
                        arcs_left_slant.append([arc_topx, arc_topy, arcbottomx,arcbottomy])
                    elif righthanded==False and slant == False:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]-50
                        arc_topy  = top_bond[1]-20
                        arcbottomx= bottom_bond[0]+50
                        arcbottomy= bottom_bond[1]
                        arcs_left.append([arc_topx, arc_topy, arcbottomx,arcbottomy])
                    elif righthanded==True and slant == True:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]
                        arc_topy  = top_bond[1]
                        arcbottomx= bottom_bond[0]-50
                        arcbottomy= bottom_bond[1]
                        arcs_right_slant.append([arc_topx, arc_topy, arcbottomx,arcbottomy])
                    elif righthanded==True and slant == False:
                        top_bond = getcoordinates[2]
                        arc_topx  = top_bond[0]+50
                        arc_topy  = top_bond[1]-20
                        arcbottomx= bottom_bond[0]-50
                        arcbottomy= bottom_bond[1]
                        arcs_right.append([arc_topx, arc_topy, arcbottomx,arcbottomy])
                    if Extra_bond==True:
                        extrabondx1=top_bond[0]
                        extrabondy1=top_bond[1]-20
                        extrabondx2=top_bond[0]
                        extrabondy2=top_bond[1]+20
                        extra_bond = [extrabondx1,extrabondy1,extrabondx2,extrabondy2]
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

                #print(All_positions_and_chains)
                #ADCs += dictionaries_to_append[10]


##Get labels and positions
            if "H[" not in keyslist[i] and "Linker[" not in keyslist[i] and "X" not in keyslist[i]:
                print("checkpoint27")
                Label_Locations = getcoordinates[3]
                text_coordinates.append([Label_Locations])
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [getcoordinates[0], righthanded,direction]
                Domain_Text.append(str(Domain_name)+mod_label)


            elif "X" in keyslist[i]:
                if mod != "Leucine":
                    ADCs.append(getcoordinates[0])
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

            if dictionary != VHa_chain and dictionary != VHb_chain and dictionary != VLa_chain and dictionary != VLb_chain and dictionary != fragment1 and dictionary != fragment2:

                if len(dictionary.get(keyslist[i])) > 1 and interaction_counter ==0:
                    if dictionary == VHa_1_test:
                        keyslistb = list(VHb_chain.keys())
                        for x in range(len(keyslistb)):
                            print(dictionary.get(keyslist[i])[1], VHb_chain.get(keyslistb[x])[0])
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


##Sort extra disulphide bridges
            if  dictionary == VHa_chain or dictionary == VHb_chain or dictionary == VLa_chain or dictionary == VLb_chain or dictionary == fragment1 or dictionary == fragment2:
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
                elif i+1 == len(dictionary) and keyslist:
                    if slant == True:
                        top_bond[1] = top_bond[1]+5
                        hinges.append([bottomx,bottomy,topx,topy])
                    elif slant == False:
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
        return(coordinates_list_heavy_a,coordinates_list_light_a,coordinates_list_heavy_b,coordinates_list_light_b,coordinates_list_heavy_c,coordinates_list_light_c,coordinates_list_heavy_d,coordinates_list_light_d,allbonds,Location_Text,text_coordinates,disulphidebridge1, disulphidebridge2 ,disulphidebridge3,disulphidebridge4,disulphidebridge5,completed_disulphidebridges,Domain_Text,Notes,Notes_positions, arcs_left,arcs_right, arcs_left_slant, arcs_right_slant,ADCs,first_interaction, hinges, linkers,names_list_heavy_a,names_list_light_a,names_list_heavy_b,names_list_light_b,names_list_heavy_c,names_list_light_c,names_list_heavy_d,names_list_light_d)

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

    if chain_count == 1:
        VHa_startx, VHa_starty = 350,100
        VHb_startx, VHb_starty = 0,0
        VLa_startx, VLa_starty = 0,0
        VLb_startx, VLb_starty = 0,0
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)

    elif  chain_count >= 3 or (chain_count == 2  and tangle_found == False and ("H[" in str(VHa_chain) and "H[" in str(VHb_chain)) or ("X[" in str(VHa_chain) and "X[" in str(VHb_chain))):
        print(str(VHa_chain))
        VHa_1_test = VHa_chain.copy()
        VLa_1_test = VLa_chain.copy()
        VHb_1_test = VHb_chain.copy()
        VLb_1_test = VLb_chain.copy()



##VHa_chain
        if "H[" in str(VHa_1_test) and "H[" in str(VHb_1_test):
            VHa_H_coordinatesx = 295
            VHa_H_coordinatesy = 280
            VHb_H_coordinatesx = 420
            VHb_H_coordinatesy = 280
        elif ("X[" in str(VHa_1_test) and "X[" in str(VHb_1_test)) and ("H[" not in str(VHa_1_test) and "H[" not in str(VHb_1_test)):
            keyslista = list(VHa_chain.keys())
            keyslistb = list(VHb_chain.keys())
            Xa = ""
            Xb = ""
            for i in range(len(keyslista)):
                if "X[" in keyslista[i]:
                    Xa = int(VHa_chain.get(keyslista[i])[0][0])
                    for i in range(len(keyslistb)):
                        if "X[" in keyslistb[i]:
                            Xb = int(VHb_chain.get(keyslistb[i])[0][0])
                            if Xa == Xb:
                                VHa_H_coordinatesx = 350
                                VHa_H_coordinatesy = 280
                                VHb_H_coordinatesx = 350
                                VHb_H_coordinatesy = 280
                            elif Xa != Xb:
                                VHa_H_coordinatesx = 295
                                VHa_H_coordinatesy = 280
                                VHb_H_coordinatesx = 420
                                VHb_H_coordinatesy = 280



        teststartx = 0
        teststarty = 0
        testHpositionVHa = renderchains(VHa_1_test,teststartx,teststarty)[25]
        print("TESTHpositive",testHpositionVHa)
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
        print("TESTHpositive", testHpositionVHb)
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




###Render Heavy chains
        print(VHa_chain)
        print(VHb_chain)
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)


###Get start positions of light chains and render
    elif chain_count == 2 :

        keyslista = list(VHa_chain.keys())
        keyslistb = list(VHb_chain.keys())
        keyslist = list(VHa_chain.keys())
        VHb_startx, VHb_starty = 400,100
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VHa_list = list(VHa_chain.keys())
        VHa_startx, VHa_starty= 200,100
        try:
            VHa_inter = VHa_chain.get(VHa_list[0])[0][1]
            print(VHa_chain.get(keyslista[0])[0][0], VHb_chain.get(keyslistb[0])[1])
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
                    VHa_startx, VHa_starty= 150,100
        except IndexError:
            if VHa_chain.get(VHa_list[1])[0][0] == VHb_chain.get(keyslistb[1])[1]:
                VHa_startx, VHa_starty= 340,100

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
                        print(VLa_startx,VLa_starty)
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
                        teststartx = 0
                        teststarty = 0
                        testHpositionVLb = renderchains(VLb_1_test,teststartx,teststarty)[25]
                        test_H_positionx = testHpositionVLb[0]
                        test_H_positiony = testHpositionVLb[1]
                        differencetest_desiredx = VLb_refx - test_H_positionx
                        differencetest_desiredy = VLb_refy - test_H_positiony
                        VLb_startx = teststartx + differencetest_desiredx
                        VLb_starty = teststarty + differencetest_desiredy
                        print(VLb_startx,VLb_starty)
                        break
    else:
        VLb_startx,VLb_starty = 0,0







    VLa_stats = renderchains(VLa_chain,VLa_startx, VLa_starty)
    VLb_stats = renderchains(VLb_chain,VLb_startx, VLb_starty)
    print("VLA_START",VLa_startx, VLa_starty)
    print("VLB_START", VLb_startx, VLb_starty)

##Two-chained Abs


    frag1_startx,frag1_starty,frag2_startx,frag2_starty,frag3_startx,frag3_starty,frag4_startx,frag4_starty = 0,0,0,0,0,0,0,0

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
        fragment_inter = fragment2.get(fragment4_list[0])[0][1]
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
    print(All_positions_and_chains)

    keyslist_all_positions_and_chains = list(All_positions_and_chains.keys())
    out_of_range = False
    how_much = 0
    for i in range(len(keyslist_all_positions_and_chains)):
        for j in range(len(All_positions_and_chains.get(keyslist_all_positions_and_chains[i])[0])):
            if isinstance(All_positions_and_chains.get(keyslist_all_positions_and_chains[i])[0][j], int) == True:
                if All_positions_and_chains.get(keyslist_all_positions_and_chains[i])[0][j] < 0:
                    out_of_range = True
                    if how_much > All_positions_and_chains.get(keyslist_all_positions_and_chains[i])[0][j]:
                        how_much = All_positions_and_chains.get(keyslist_all_positions_and_chains[i])[0][j]

    if out_of_range == True:
        new_start = how_much+5




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

    if out_of_range == True:
        new_start = how_much-10
        coordinates_to_change = [Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Bonds,Hinges,Linkers,ADCs,disulphide_bridges, Label_spot]
        for i in range(len(coordinates_to_change)):
            for j in range(len(coordinates_to_change[i])):
                for k in range(len(coordinates_to_change[i][j])):
                    print(coordinates_to_change[i][j][k])
                    if isinstance(coordinates_to_change[i][j][k], int) == True or isinstance(coordinates_to_change[i][j][k], float) == True:
                        if  k %2 != 0:
                            coordinates_to_change[i][j][k] -= new_start
                    else:
                        for l in range(len(coordinates_to_change[i][j][k])):
                            if isinstance(coordinates_to_change[i][j][k][l], int) == True or isinstance(coordinates_to_change[i][j][k][l], float) == True:
                                if l %2 != 0:
                                    coordinates_to_change[i][j][k][l] -= new_start

    return(Bonds,disulphide_bridges,Hinges,Linkers,Heavy_Domains_a,names_list_heavy_a,Light_Domains_a,names_list_light_a,Heavy_Domains_b,names_list_heavy_b,Light_Domains_b,names_list_light_b,Heavy_Domains_c,names_list_heavy_c,Light_Domains_c,names_list_light_c,Heavy_Domains_d,names_list_heavy_d,Light_Domains_d,names_list_light_d,Label_Text,Label_spot,Domain_Text,Notes,Notes_positions,arcs_left,arcs_right,arcs_left_slant,arcs_right_slant,ADCs)

def render(chains_list,canvas,text_to_image):
    if text_to_image == True:
        canvas.delete("all")
        global canvas_polygons
        canvas_polygons = {}
        global canvas_labels
        canvas_labels= {}
        global Label_lock
        global custom_labels
        custom_labels = {}

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

        print(Heavy_Domains_a, names_Heavy_a)
        print(Light_Domains_a, names_Light_a)

    #disulphide_bridge
        if disulphide_bridges != []:
            for i in range(len(disulphide_bridges)):
                domain = canvas.create_line(disulphide_bridges[i], fill='#FF4040', width = 2,tags="disulphide")
                canvas_polygons[domain] = [disulphide_bridges[i], "-"]

    #Bonds
        for i in range(len(Bonds)):
            domain = canvas.create_line(Bonds[i], fill='#000000', width = 2,tags="bonds")
            canvas_polygons[domain] = [Bonds[i], "-"]
        if arcs_left!=[]:
            for i in range(len(arcs_left)):
                domain = canvas.create_arc(arcs_left[i], start=90, extent=180, style=tk.ARC, fill='#000000', width = 2,tags="bonds")
                canvas_polygons[domain] = [arcs_left[i], "-"]
        if arcs_right!=[]:
            for i in range(len(arcs_right)):
                domain = canvas.create_arc(arcs_right[i], start=270, extent=180, style=tk.ARC, fill='#000000', width = 2,tags="bonds")
                canvas_polygons[domain] = [arcs_right[i], "-"]
        if arcs_left_slant != []:
            for i in range(len(arcs_left_slant)):
                domain = canvas.create_arc(arcs_left_slant[i], start=150, extent=120, style=tk.ARC,width=2,tags="bonds")
                canvas_polygons[domain] = [arcs_left_slant[i], "-"]
        if arcs_right_slant != []:
            for i in range(len(arcs_right_slant)):
                domain = canvas.create_arc(arcs_right_slant[i], start=270, extent=120, style=tk.ARC,width=2,tags="bonds")
                canvas_polygons[domain] = [arcs_right_slant[i], "-"]
        for i in range(len(Linkers)):
            domain = canvas.create_line(Linkers[i], fill='#000000', width = 2,tags="bonds")
            canvas_polygons[domain] = [Linkers[i], "-L-"]
        for i in range(len(Hinges)):
            domain = canvas.create_line(Hinges[i], fill='#000000', width = 2,tags="bonds")
            canvas_polygons[domain] = [Hinges[i], "-H-"]

    #A domains
        for i in range(len(Heavy_Domains_a)):
            domain = canvas.create_polygon(Heavy_Domains_a[i], outline='#000000',fill='#007ECB', width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_a[i], names_Heavy_a[i]]
        for i in range(len(Light_Domains_a)):
            domain = canvas.create_polygon(Light_Domains_a[i], outline='#000000',fill='#73CAFF', width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_a[i], names_Light_a[i]]
    #B domains
        for i in range(len(Heavy_Domains_b)):
            domain = canvas.create_polygon(Heavy_Domains_b[i], outline='#000000',fill='#FF43EE', width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_b[i], names_Heavy_b[i]]
        for i in range(len(Light_Domains_b)):
            domain = canvas.create_polygon(Light_Domains_b[i], outline='#000000',fill='#F9D3F5', width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_b[i], names_Light_b[i]]
    #C domains
        for i in range(len(Heavy_Domains_c)):
            domain = canvas.create_polygon(Heavy_Domains_c[i], outline='#000000',fill='#0BD05A', width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_c[i], names_Heavy_c[i]]
        for i in range(len(Light_Domains_c)):
            domain = canvas.create_polygon(Light_Domains_c[i], outline='#000000',fill='#B9FAD3', width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_c[i], names_Light_c[i]]
    #D domains
        for i in range(len(Heavy_Domains_d)):
            domain = canvas.create_polygon(Heavy_Domains_d[i], outline='#000000',fill='#D9DE4A', width=2,tags="domain")
            canvas_polygons[domain] = [Heavy_Domains_d[i], names_Heavy_d[i]]
        for i in range(len(Light_Domains_d)):
            domain = canvas.create_polygon(Light_Domains_d[i], outline='#000000',fill='#E2F562', width=2,tags="domain")
            canvas_polygons[domain] = [Light_Domains_d[i], names_Light_d[i]]
    #ADCs
        if ADCs != []:
            for i in range(len(ADCs)):
                domain = canvas.create_polygon(ADCs[i], outline='#000000',fill='#68C1C1', width=2,tags="domain")
                canvas_polygons[domain] = [ADCs[i],  "X"]
    #Labels
        if Label_lock == True:
            for i in range(len(Label_positions)):
                x = Label_positions[i][0][0]
                y = Label_positions[i][0][1]
                print(x,y)
                Domain_Text[i] = re.sub("sd","",Domain_Text[i])
                label  = canvas.create_text(x,y, text=Domain_Text[i],tags = "label")
                canvas_labels[label] = [[x,y], Domain_Text[i]]


        if Notes != []:
            print(Notes)
            notes_set = set(Notes)
            setlist = list(notes_set)
            for i in range(len(setlist)):
                note = canvas.create_text(Note_positions[i],text=setlist[i],tags = "label")
                custom_labels[note] = [Note_positions[i],setlist[i]]


        print(canvas_polygons)
        print(canvas_labels)



def sequence_render_pipeline(canvas):
    '''
    Gets sequence and automatically tidies it
    '''
    global Bond_lock
    global Delete_lock
    Delete_lock = False
    Bond_lock = ""
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
    Delete_lock = False
    Bond_lock = ""
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
    elif Domain_Primer_Lock == domain_name:
        Domain_Primer_Lock = ""
        Bond_lock = ""
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        lower_canvas.config(cursor = "arrow")
        status_label.config(text="")
        Domain_Primer = ()
    print(Domain_Primer)
############################################

def domain_button(canvas,startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H,domain_name,Light,Heavy):
    '''
    Draw domains onto canvas with domainmaker
    '''
    global Bond_lock
    global Delete_lock
    global Label_lock
    Delete_lock = False
    Bond_lock = ""
    lower_canvas.bind("<Button-1>", mm.place_domain)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "arrow")
    status_label.config(text="")
    domaincoordinates = domainmaker(startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H)
    if "a" in str(domain_name):
        heavy_colour, light_colour = '#007ECB', '#73CAFF'
    elif "b" in str(domain_name):
        heavy_colour, light_colour = '#FF43EE', '#F9D3F5'
    elif "c" in str(domain_name):
        heavy_colour, light_colour = '#0BD05A', '#B9FAD3'
    elif "d" in str(domain_name):
        heavy_colour, light_colour = '#D9DE4A', '#E2F562'
    elif "X" in str(domain_name):
        heavy_colour, light_colour = '#68C1C1','#68C1C1'
    else:
        heavy_colour, light_colour = '#5C5C5C','#B0B0B0'
    if Heavy == True:
        domain = lower_canvas.create_polygon(domaincoordinates[0], outline='#000000',fill=heavy_colour, width=2, tags="domain")
    elif Light == True:
        domain = lower_canvas.create_polygon(domaincoordinates[0], outline='#000000',fill=light_colour, width=2, tags="domain")
    canvas_polygons[domain] = [domaincoordinates[0], domain_name]
    if Label_lock == True:
        domain_name = re.sub("\.|@|>","",domain_name)
        label  = lower_canvas.create_text(domaincoordinates[3], text = str(domain_name), tags = "label")
        canvas_labels[label] = [domaincoordinates[3], domain_name]
    global domain_mod
    domain_mod = ""
    global domain_direction
    domain_direction = "constant"

############################################
def save_as_png(canvas):
    fileName = "AbYdraw_export"
    # save postscipt image
    canvas.postscript(file = fileName + '.eps')
    # use PIL to convert to PNG
    img = Image.open(fileName + '.eps')
    #img.save(fileName + '.png', 'png')
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
    Delete_lock = False
    Bond_lock = ""
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    lower_canvas.config(cursor = "arrow")
    status_label.config(text="")
    output = open("Template_File_Export.txt","w")
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
        LVJC_Germline_Strings = []
        L_chain_Ranges_Strings = []
        H_CDR_Strings = []
        L_CDR_Strings =[]
        all = [Type,Chain_Class_Strings,NGlycos_Strings,Heavy_CysPositions,Heavy_Disulphides_Intra,Light_CysPositions,Light_Disulphides_Intra,DisulphidesInter,HVJC_Germline_Strings,H_chain_Ranges_Strings,LVJC_Germline_Strings,L_chain_Ranges_Strings,H_CDR_Strings,L_CDR_Strings, [], []]
        index = str("["+str(indexing)+"]")
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
            listier = re.sub("\[[0-9]\]|\[|\]|'|,","",str(listing))
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
        if "VH" and "VL" in str(listing):
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
    for i in Template_File:
        output.write(str(i)+"\n")





############################################
def sequence_pipeline(canvas):
    '''
    Take drawing from cavnas and generate expression to be displayed on entry box
    '''
    lower_canvas.config(cursor = "arrow")
    Bond_lock = ""
    Delete_lock = False
    lower_canvas.bind("<Button-1>", mm.select)
    lower_canvas.bind("<B1-Motion>", mm.drag)
    lower_canvas.bind("<ButtonRelease-1>", mm.release)
    polygons_keyslist = list(canvas_polygons.keys())
    custom_keyslist = list(custom_labels.keys())
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
    print(disulphides_list)
    disulphides_dict = {}
    for i in range(len(disulphides_list)):
        for j in range(len(polygons_keyslist)):
            if disulphides_list[i] == polygons_keyslist[j]:
                disulphides_dict[j] = canvas_polygons.get(polygons_keyslist[j])
    custom_list = canvas.find_withtag("custom_labels")
    custom_dict = {}
    for i in range(len(custom_list)):
        for j in range(len(custom_keyslist)):
            if custom_list[i] == custom_keyslist[j]:
                custom_dict[j] = custom_labels.get(custom_keyslist[j])

    print("DOMAINs", domains_dict)
    print("BONDS", bonds_dict)
    print("Disulphides", disulphides_dict)
    print("Labels", custom_dict)
    chains=[]
    current_chain_str = []
    current_chain_coords_lists = []
    #directions = []
    domains_keyslist = list(domains_dict.keys())
    bonds_keyslist = list(bonds_dict.keys())
    disulphides_keyslist = list(disulphides_dict.keys())
    labels_keyslist = list(custom_dict.keys())
    print(labels_keyslist)

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
                ##find domain that bond is attached to
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
                            if domainx1 < bondx1 < domainx2 and domainy1 < bondy1 < domainy2:
                                connection_found = True
                                full_chain.append(bonds_keyslist[n])
                                string.append(bonds_dict.get(bonds_keyslist[n])[1])
                                #directions.append(bonds_dict.get(bonds_keyslist[n])[2])
                                break

                if connection_found == True:
                    continuation_found = True
                else:
                    continuation_found = False


        #full_directions.append(directions)
        full_chains.append(full_chain)
        strings.append(string)

##number chains
    counter = 1
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            if str(strings[i][j]) != "-" and  str(strings[i][j]) != "-L-":
                if "-" not in strings[i][j]:
                    strings[i][j] += str("("+str(counter)+")")
                elif "-H-" in strings[i][j]:
                    strings[i][j] = str("-H("+str(counter)+")-")
                elif "-X-" in strings[i][j]:
                    strings[i][j] = str("-X("+str(counter)+")-")
                counter += 1


    print(strings)
##Pair chains based on closeness
    paired = []
    for i in range(len(strings)):
        for j in range(len(strings[i])):
            if ":" not in str(strings[i][j]) and "-" not in str(strings[i][j]) and "sd" not in str(strings[i][j]):
                number =  re.findall("\((.*?)\)", str(strings[i][j]))
                number =  int(re.sub("\[|\'|\]","", str(number)))

                if number not in paired:
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
                                                print(d1x1, g[0],d1x2,"HMM", d1y1, g[1],d1y2)
                                                if d1x1 < g[0] < d1x2 and d1y1 < g[1] < d1y2:
                                                    print("YES ALRIGHT")
                                                    if ("VH" in str(strings[i][j]) and "VL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("VL" in str(strings[i][j]) and "VH" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH1" in str(strings[i][j]) and "CH1" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CL" in str(strings[i][j]) and "CH1" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH1" in str(strings[i][j]) and "CL" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH2" in str(strings[i][j]) and "CH2" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH3" in str(strings[i][j]) and "CH3" in str(domains_dict.get(domains_keyslist[f])[1])) or ("CH4" in str(strings[i][j]) and "CH4" in str(domains_dict.get(domains_keyslist[f])[1])) or ("-H-" == str(strings[i][j]) and "-H-" == str(domains_dict.get(domains_keyslist[f])[1])) or ("X" in str(strings[i][j]) and "X" in str(domains_dict.get(domains_keyslist[f])[1])):
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
                                            print("H1", d1x1,d1x2,d1y1,d1y2," H2 ",d2x1,d2x2,d2y1,d2y2)
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



    for i in range(len(full_chains)):
        for j in range(len(full_chains[i])):
            #print(domains_dict.get(full_chains[i][j])[0])
            if "-" not in strings[i][j]:
                coordinates = domains_dict.get(full_chains[i][j])[0]
                min_max = get_min_max_coordinates(coordinates)
                d1x1 = min_max[0]
                d1x2 = min_max[1]
                d1y1 = min_max[2]
                d1y2 = min_max[3]
                print(d1x1,d1x2,d1y1,d1y2)
                for k in range(len(labels_keyslist)):
                    labelx = custom_dict.get(labels_keyslist[k])[0][0]
                    labely = custom_dict.get(labels_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if (d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2):
                        note = custom_dict.get(labels_keyslist[k])[1]
                        noting = str("[TYPE:"+note+"]")
                        strings[i][j] += noting
            elif "-H" in strings[i][j]:
                coordinates = bonds_dict.get(full_chains[i][j])[0]
                min_max = get_min_max_coordinates(coordinates)
                d1x1 = min_max[0]-50
                d1x2 = min_max[1]+50
                d1y1 = min_max[2]
                d1y2 = min_max[3]
                print(d1x1,d1x2,d1y1,d1y2)
                for k in range(len(labels_keyslist)):
                    labelx = custom_dict.get(labels_keyslist[k])[0][0]
                    labely = custom_dict.get(labels_keyslist[k])[0][1]
                    #print(labelx,labely)
                    if (d1x1 <= labelx <= d1x2) and (d1y1 <= labely <= d1y2):
                        domain = strings[i][j].split("-")[1]
                        note = custom_dict.get(labels_keyslist[k])[1]
                        noting = str("-"+domain+"[TYPE:"+note+"]-")
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
    final_string = re.sub("\-L\-\|","|",str(final_string))
    final_string = re.sub("\-\|","|",str(final_string))


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
    global custom_labels
    canvas_polygons = {}
    canvas_labels = {}
    temp_label = {}
    custom_labels = {}

def delete_button(canvas):
    '''
    Mouse click deletes object from canvas
    '''
    global Delete_lock
    global Bond_lock
    Bond_lock = ""
    if Delete_lock == False:
        Delete_lock = True
        lower_canvas.config(cursor = "arrow")
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>",mm.delete)
        status_label.config(text="Delete lock on")

    elif Delete_lock == True:
        Delete_lock = False
        lower_canvas.unbind("<Button-1>")
        lower_canvas.bind("<Button-1>", mm.select)
        lower_canvas.bind("<B1-Motion>", mm.drag)
        lower_canvas.bind("<ButtonRelease-1>", mm.release)
        status_label.config(text="")

def domain_type_button(letter):
    domain_letter = str(letter)
    global domain_type
    domain_type = letter
def domain_mod_button(letter):
    domain_letter = str(letter)
    global domain_mod
    domain_mod = str(letter)
    if letter == ">" or letter == "@":
        global domain_direction
        domain_direction = "innie"

def bond_drag_button(canvas,name,buttonpress):
    global Bond_lock
    global Delete_lock
    Delete_lock = False
    if Bond_lock != buttonpress:
        Bond_lock = buttonpress
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.start_bond)
        lower_canvas.bind("<B1-Motion>", mm.drag_bond)
        if name == "-":
            lower_canvas.bind("<ButtonRelease-1>", mm.release_bond)
            status_label.config(text="Bond lock on")
        elif name == "-H-":
            lower_canvas.bind("<ButtonRelease-1>", mm.release_Hinge_bond)
            status_label.config(text="Hinge lock on")
        elif name == "-L-":
            lower_canvas.bind("<ButtonRelease-1>", mm.release_Linker_bond)
            status_label.config(text="Linker lock on")
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
    Delete_lock = False
    if Bond_lock != buttonpress:
        Bond_lock = buttonpress
        lower_canvas.unbind("<Button-1>")
        lower_canvas.unbind("<B1-Motion>")
        lower_canvas.unbind("<ButtonRelease-1>")
        lower_canvas.bind("<Button-1>", mm.start_bond)
        lower_canvas.bind("<B1-Motion>", mm.drag_disulphide_bond)
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

def CustomLabelButton_function(canvas):
        global CustomLabelLock
        if CustomLabelLock == False:
            CustomLabelLock = True
            entry=CustomLabelEntry.get("1.0","end-1c")
            lower_canvas.bind("<Button-1>", mm.place_custom_label)
            lower_canvas.bind("<ButtonRelease-1>", mm.place_domain_release)
            status_label.config(text=entry)
            lower_canvas.config(cursor = "plus")

        elif CustomLabelLock == True:
            CustomLabelLock = False
            lower_canvas.unbind("<Button-1>")
            lower_canvas.unbind("<B1-Motion>")
            lower_canvas.unbind("<ButtonRelease-1>")
            lower_canvas.bind("<Button-1>", mm.select)
            lower_canvas.bind("<B1-Motion>", mm.drag)
            lower_canvas.bind("<ButtonRelease-1>", mm.release)
            status_label.config(text="")
            lower_canvas.config(cursor = "arrow")


def items_selected(e):
    '''
    Render item selected in library
    '''
    textBox.delete("1.0","end")
    global Bond_lock
    global Delete_lock
    global antibodyformats
    global formats_keyslist
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
        global custom_labels
        label_keyslist = list(canvas_labels.keys())
        custom_keyslist = list(custom_labels.keys())
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

        elif self.item in custom_keyslist:
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

                    if x1<= labelx <=x2 and y1 <= labely <= y2:
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
        custom_keyslist = list(custom_labels.keys())
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

        elif self.item in custom_keyslist:
            x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
            y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
            diffx = x2-x1
            diffy = y2-y1
            coordinates    = custom_labels.get(self.item)[0]
            name           = custom_labels.get(self.item)[1]
            new_coordinates= []
            for i in range(len(coordinates)):
                if i%2 ==0:
                    new_coordinates.append((coordinates[i]+diffx))
                elif i%2!=0:
                    new_coordinates.append((coordinates[i]+diffy))
            custom_labels[self.item]=[new_coordinates, name]
            print(custom_labels.get(self.item))
        print(canvas_polygons)
        print(canvas_labels)
        startcoordinates = []
        newcoordinates = []
    ####Delete selected item on canvas###
    def delete(self,event):
        global canvas_polygons
        global canvas_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc)[0]        # ID for closest
        polygons_keyslist = list(canvas_polygons.keys())
        labels_keyslist = list(canvas_labels.keys())
        if self.item in labels_keyslist:
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
        min_max = get_min_max_coordinates(domain_coordinates)
        x1 = min_max[0]
        x2 = min_max[1]
        y1 = min_max[2]
        y2 = min_max[3]
        for i in range(len(labels_keyslist)):
            labelx = canvas_labels.get(labels_keyslist[i])[0][0]
            labely = canvas_labels.get(labels_keyslist[i])[0][1]
            if x1< labelx <x2 and y1 < labely < y2:
                lower_canvas.delete(labels_keyslist[i])
                del canvas_labels[labels_keyslist[i]]
        lower_canvas.delete(self.item)
        del canvas_polygons[self.item]
    ###Click item to reverse orientation###
    def change_orientation(self,event):
        self.startcoordinates = []
        self.newcoordinates =[]
        widget = event.widget                       # Get handle to canvas
        # Convert screen coordinates to canvas coordinates
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        self.item = widget.find_closest(xc, yc,halo = 5, start="domain")[0]
        global canvas_labels
        global temp_label
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

            if "a" in str(domain_name):
                heavy_colour, light_colour = '#007ECB', '#73CAFF'
            elif "b" in str(domain_name):
                heavy_colour, light_colour = '#FF43EE', '#F9D3F5'
            elif "c" in str(domain_name):
                heavy_colour, light_colour = '#0BD05A', '#B9FAD3'
            elif "d" in str(domain_name):
                heavy_colour, light_colour = '#D9DE4A', '#E2F562'
            else:
                heavy_colour, light_colour = '#5C5C5C','#B0B0B0'
            if "VH" in domain_name or "CH" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=heavy_colour, width=2, tags="domain")
            elif "VL" in domain_name or "CL" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill=light_colour, width=2, tags="domain")
            elif "X" in domain_name:
                domain = lower_canvas.create_polygon(new_coordinates, outline='#000000',fill='#68C1C1', width=2, tags="domain")

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
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        startx = xc
        starty = yc-40
        Domain_Primer[5] = re.sub("\+|\-","",Domain_Primer[5])
        domaincoordinates = domainmaker(startx,starty,Domain_Primer[0],Domain_Primer[1],Domain_Primer[2],Domain_Primer[3],Domain_Primer[4],Domain_Primer[5],Domain_Primer[6],Domain_Primer[7])
        domain_name = Domain_Primer[8]
        if "a" in str(domain_name):
            heavy_colour, light_colour = '#007ECB', '#73CAFF'
        elif "b" in str(domain_name):
            heavy_colour, light_colour = '#FF43EE', '#F9D3F5'
        elif "c" in str(domain_name):
            heavy_colour, light_colour = '#0BD05A', '#B9FAD3'
        elif "d" in str(domain_name):
            heavy_colour, light_colour = '#D9DE4A', '#E2F562'
        elif "X" in str(domain_name):
            heavy_colour, light_colour = '#68C1C1','#68C1C1'
        else:
            heavy_colour, light_colour = '#5C5C5C','#B0B0B0'
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
        domain_mod = ""
        global domain_direction
        domain_direction = "constant"

    def place_domain_release(self,event):
        lower_canvas.config(cursor = "plus")

    def place_custom_label(self,event):
        global custom_labels
        widget = event.widget
        xc = widget.canvasx(event.x); yc = widget.canvasx(event.y)
        entry=CustomLabelEntry.get("1.0","end-1c")
        if entry != "":
            label = lower_canvas.create_text(xc,yc, text = entry, tags = "custom_labels")
            custom_labels[label] = [[xc,yc], entry]
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
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#000000', width=2, tags="bonds")
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates   = []
    def release_Hinge_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-H-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#000000', width=2, tags="bonds")
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordiantes   = []
    def release_Linker_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-L-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#000000', width=2, tags="bonds")
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
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#FF4040', width=2,tags = "draggable_line")
        self.newcoordinates = [x2,y2]
    def release_Disulphide_bond(self,event):
        lower_canvas.delete("draggable_line")
        x1,x2 = self.startcoordinates[0], self.newcoordinates[0]
        y1,y2 = self.startcoordinates[1], self.newcoordinates[1]
        name = "-"
        domain = lower_canvas.create_line(x1,y1,x2,y2, fill='#FF4040', width=2, tags="disulphide")
        canvas_polygons[domain] = [[x1,y1,x2,y2],name]
        self.startcoordinates = []
        self.newcoordinates =[]










################Domain Drawer######################
def domainmaker(startx,starty,righthanded,slant,V,direction,X,mod,interaction,previous_H):
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
domain_type = "a"
domain_mod  = ""
domain_direction = "constant"
a_button= tk.Button(frame2,text="a",bg = "grey", command =lambda: domain_type_button("a"))
a_button.place(relx = 0.41, rely = 0.21, relheight = 0.1, relwidth=0.1)
b_button= tk.Button(frame2,text="b",bg = "grey", command =lambda: domain_type_button("b"))
b_button.place(relx = 0.51, rely = 0.21, relheight = 0.1, relwidth=0.1)
c_button= tk.Button(frame2,text="c",bg = "grey", command =lambda: domain_type_button("c"))
c_button.place(relx = 0.41, rely = 0.31, relheight = 0.1, relwidth=0.1)
d_button= tk.Button(frame2,text="d",bg = "grey", command =lambda: domain_type_button("d"))
d_button.place(relx = 0.51, rely = 0.31, relheight = 0.1, relwidth=0.1)
KIH_knob= tk.Button(frame2,text="Knob",bg = "grey", command =lambda: domain_mod_button(">"))
KIH_knob.place(relx = 0.41, rely = 0.41, relheight = 0.1, relwidth=0.1)
KIH_hole= tk.Button(frame2,text="Hole",bg = "grey", command =lambda: domain_mod_button("@"))
KIH_hole.place(relx = 0.51, rely = 0.41, relheight = 0.1, relwidth=0.1)
Positive_charge= tk.Button(frame2,text="+",bg = "grey", command =lambda: domain_mod_button("+"))
Positive_charge.place(relx = 0.41, rely = 0.51, relheight = 0.1, relwidth=0.1)
Negative_charge= tk.Button(frame2,text="-",bg = "grey", command =lambda: domain_mod_button("-"))
Negative_charge.place(relx = 0.51, rely = 0.51, relheight = 0.1, relwidth=0.1)
Zip_button=tk.Button(frame2,text="Zip",bg = "grey", command =lambda: domain_mod_button("Leucine"))
Zip_button.place(relx = 0.41, rely = 0.61, relheight = 0.1, relwidth=0.1)
sdFV_button = tk.Button(frame2,text="sdFV",bg = "grey", command =lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"Single_Fv_Chain",False,domain_mod,"","",str("sdVH"+domain_mod+"."+domain_type),False,True))
sdFV_button.place(relx = 0.41, rely = 0.01, relheight = 0.2, relwidth=0.2)
###Insert bonds buttons ###
##Col1
InsertVHDomainButton= tk.Button(frame2,text="VH",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"innie",False,domain_mod,"","",str("VH"+domain_mod+"."+domain_type),False,True))
InsertVHDomainButton.place(relx = 0.01, rely = 0.01, relheight = 0.2, relwidth=0.2)
InsertCH1DomainButton= tk.Button(frame2,text="CH1",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH1"+domain_mod),False,True))
InsertCH1DomainButton.place(relx = 0.01, rely = 0.21, relheight = 0.2, relwidth=0.2)
InsertCH2DomainButton= tk.Button(frame2,text="CH2",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH2"+domain_mod),False,True))
InsertCH2DomainButton.place(relx = 0.01, rely = 0.41, relheight = 0.2, relwidth=0.2)
InsertCH3DomainButton= tk.Button(frame2,text="CH3",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH3"+domain_mod),False,True))
InsertCH3DomainButton.place(relx = 0.01, rely = 0.61, relheight = 0.2, relwidth=0.2)
##Col2
InsertVLDomainButton= tk.Button(frame2,text="VL",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,True,"outie",False,domain_mod,"","",str("VL"+domain_mod+"."+domain_type),True,False))
InsertVLDomainButton.place(relx = 0.21, rely = 0.01, relheight = 0.2, relwidth=0.2)
InsertCLDomainButton= tk.Button(frame2,text="CL",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CL"+domain_mod),True,False))
InsertCLDomainButton.place(relx = 0.21, rely = 0.21, relheight = 0.2, relwidth=0.2)
InsertCH4DomainButton= tk.Button(frame2,text="CH4",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,False,domain_mod,"","",str("CH4"+domain_mod),False,True))
InsertCH4DomainButton.place(relx = 0.21, rely = 0.41, relheight = 0.2, relwidth=0.2)
InsertXDomainButton= tk.Button(frame2,text="Other",bg = "grey", command=lambda: prime_domain_button(lower_canvas, 400,100,True,False,False,domain_direction,True,domain_mod,"","",str("X"+domain_mod),True,False))
InsertXDomainButton.place(relx = 0.21, rely = 0.61, relheight = 0.2, relwidth=0.2)
###Drag and pull bonds###
Bond_lock = ""
InsertBondButton= tk.Button(frame2,text="Connect",bg = "grey", command=lambda: bond_drag_button(lower_canvas,"-","bond"))
InsertBondButton.place(relx = 0.61, rely = 0.01, relheight = 0.2, relwidth=0.2)
InsertDBondButton= tk.Button(frame2,text="Disulphide",bg = "grey", command=lambda: disulphide_bond_button(lower_canvas,"-","disulphide"))
InsertDBondButton.place(relx = 0.61, rely = 0.21, relheight = 0.2, relwidth=0.2)
InsertLHingeButton= tk.Button(frame2,text="Hinge",bg = "grey", command=lambda: bond_drag_button(lower_canvas,"-H-","hinge"))
InsertLHingeButton.place(relx = 0.61, rely = 0.41, relheight = 0.2, relwidth=0.2)
InsertLinkerButton= tk.Button(frame2,text="Linker",bg = "grey", command=lambda: bond_drag_button(lower_canvas,"-L-","linker"))
InsertLinkerButton.place(relx = 0.61, rely = 0.61, relheight = 0.2, relwidth=0.2)

###Delete/clear###
Delete_lock = False
Label_lock = True
CustomLabelLock = False
InsertDelAllButton = tk.Button(frame2,text="Clear All",bg="grey", command=lambda: delete_all_button(lower_canvas))
InsertDelAllButton.place(relx = 0.81,rely = 0.01, relheight = 0.2, relwidth= 0.18)
InsertDelClickButton = tk.Button(frame2,text="Delete", bg="grey", command=lambda: delete_button(lower_canvas))
InsertDelClickButton.place(relx = 0.81,rely = 0.21, relheight = 0.2, relwidth= 0.18)
Labels_buttons = tk.Button(frame2,text="Labels", bg="grey", command=lambda: labels_button(lower_canvas))
Labels_buttons.place(relx = 0.81,rely = 0.41, relheight = 0.2, relwidth= 0.18)
CustomLabelButton = tk.Button(frame2,text="Custom Label", bg="grey", command=lambda: CustomLabelButton_function(lower_canvas))
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
"Knobs in Holes":"VH.b(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3@(5:12) | VL.b(6:1)-CL(7:2){1} | VH.a(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3>(12:5) | VL.a(13:8)-CL(14:9){1}",
"IgG-H-scFV":"VH.a(1:8)-CH1(2:9){1}-H(3:12){2}-CH2(4:13)-CH3(5:14)-L-VH.b(6:7)-L-VL.b(7:6)|VL.a(8:1)-CL(9:2){1}|VH.a(10:17)-CH1(11:18){1}-H(12:3){2}-CH2(13:4)-CH3(14:5)-L-VH.b(15:16)-L-VL.b(16:15)|VL.a(17:10)-CL(18:11){1}",
"IgG-L-scFV":"VH.b(1:6) -CH1(2:7){1} -H(3:12){2} -CH2(4:13) -CH3(5:14) |VL.b(6:1) -CL(7:2){1} -L -VH.b(8:9) -L -VL.b(9:8) |VH.a(10:15) -CH1(11:16){1} -H(12:3){2} -CH2(13:4) -CH3(14:5) | VL.a(15:10) -CL(16:11){1} -L -VH.a(17:18) -L -VL.a(18:17)",
"scFV-H-scFV": "VL.b(1:2)-L-VH.b(2:1)-H(3:8){2}-VH.c(4:5)-L-VL.c(5:4)|VL.a(6:7)-L-VH.a(7:6)-H(8:3){2}-VH.d(9:10)-L-VL.d(10:9)",
"F(ab)2":"VH.b(1:4) - CH1(2:5){1} -H(3:8){2} | VL.b(4:1) - CL(5:2){1} | VH.a(6:9) - CH1(7:10){1} - H(8:3){2} | VL.a(9:6) - CL(10:7){1}",
"scFV":"VH.a(1:2)-L-VL.a(2:1)",
"scFV4":"VH.b(4:5)-L-VL.b(5:4)-CH1(6:3){1}-H(7:16){2}-CH2(8:17)-CH3(9:18)|  VL.b(1:2)-L-VH.b(2:1)-CL(3:6){1} | VH.a(13:14)-L-VL.a(14:13)-CH1(15:12){1}-H(16:7){2}-CH2(17:8)-CH3(18:9) | VL.a(10:11)-L-VH.a(11:10)-CL(12:15){1}",
"sdFV4":"VH.b(1) -CH1(2:7){1} -H(3:10){2}-CH2(4:11) -CH3@(5:12) | VH.b(6)  -CL(7:2){1} | VH.a(8) -CH1(9:14){1}-H(10:3){2}-CH2(11:4) -CH3>(12:5) | VH.a(13) -CL(14:9){1}",
"Nanobody":"VH.a(1)",
"BiTE":"VH.a(1:2) -L -VL.a(2:1) -L -VH.b(3:4) -L -VL.b(4:3)",
"HSAbody":"VL.a(1:2)-L-VH.a(2:1)-L-X(3)[TYPE:FUSION][NOTE:human serum albumin]-L-VH.b(4:5)-L-VL.b(5:4)",
"Cov-X-body":"X(1)[TYPE: pharmacophore peptide heterodimer]-VH.b(2:7)- CH1(3:8){1}-H(4:12){2}-CH2(5:11)-CH3(6:12) | VL.b(7:2)-CL(8:3){1} | X(9)[TYPE: pharmacophore peptide heterodimer]-VH.a(10:15)-CH1(11:16){1}-H(12:4){2}- CH2(13:5)-CH3 (14:6) | VL.a(15:10)-CL(16:11){1}",
"Diabody":"VH.a(1:3) -L -VH.b(2:4) | VL.a(3:1) -L -VL.b(4:2)",
"Miniantibody":"VH.b(1:2)-L-VL.b(2:1)-H*(3:7){1}[MOD: engineered disulphide bond]-X(4)[TYPE:LEUCINE] | VH.a(5:6)-L-VL.a(6:5)-H*(7:3){1}[MOD: engineered disulphide bond]-X(8)[TYPE:LEUCINE]",
"scDiabody":"VH.a(1:4) -L -VL.a(2:3) -L -VH.b(3:2) -L -VL.b(4:1)",
"scDiabody-CH3":"VL.b(1:4){1} -L -VL.b(2:3) -L -VH.b(3:2) -L -VH.b(4:1){1}-H(5:12){2}-CH2(6:13)-CH3(7:14)| VL.a(8:11){1}-L-VL.a(9:10)-L-VH.a(10:9)-L-VH.a(11:8){1}-H(12:5){2}-CH2(13:6)-CH3(14:7)",
"DART":"VL.a(1:5) -L -VH.b(2:4) -H*(3:6){3}[MOD: engineered disulphide bond] | VL.b(4:2) -L -VH.a(5:1) -H*(6:3){3}[MOD: engineered disulphide bond]",
"Intrabody":"VL.b(1:2) -L -VH.b(2:1)-H(3:10){2}- CH1(4:11) - CH2(5:12) -L - VH.a(6:7) -L -VL.a(7:6)| VL.a(8:9) -L -VH.a(9:8) -H(10:3){2} -CH1(11:4) - CH2(12:5) -L -VL.b(13:14) -L -VH.b(14:13)",
"Fv-Fc":"VH.b(1:5){1}-H(2:7){2}-CH1(3:8)-CH2(4:9)| VL.b(5:1){1}|VH.a(6:10){1}-H(7:2){2}-CH1(8:3)-CH2(9:4)|VL.a(10:6){1}",
"Triplebody":"VH.a(1:5)-CH1(2:6){1} -L- VL.b(3:4) -L -VH.b(4:3) | VL.a(5:1) -CL(6:2){1} -L - VL.c(7:8) -L -VH.c(8:7)",
"scTriplebody":"VH.a(1:5)-CH1(2:6){2}-L-VH.b(3:4){1}-L-VL.b(4:3){1}-L-VL.a(5:1)-CH2(6:2){2}-L-VH.c(7:8){1}-L-VL.c(8:7){1}",
"TriBiMinibody":"VH.b(1:2) -L -VL.b(2:1) -H(3:9){2} -CH3@(4:10){2} -L -VH.c(5:6) -L - VL.c(6:5) | VH.a(7:8) -L -VL.a(8:7) -H(9:3){2}-CH3>(10:4)",
"LUZ-Y":"VL.b(1:3)-CL(2:4){1}-L -VH.b(3:1)-CH1(4:2){1}-H(5:13){2}-CH2(6:14) -CH3(7:15) -X(8)[TYPE: LEUCINE] | VL.a(9:11)-CL(10:12){1} -L -VH.a(11:9) -CH1(12:10){1}-H(13:5){2}-CH2(14:6)-CH3(15:7)-X(16)[TYPE: LEUCINE]",
"Dock and Lock":"VH.b(1:6)-CH1(2:7){1}-L-X(3:3)[TYPE:FUSION]-L-CH2(4:14){1}-VH.c(5:15) | VL.b(6:1)-CL(7:2){1} | VH.a(8:10)-CH1(9:13){1}-L-X(3:3)[TYPE:FUSION]-L-CH2(10:16){1}-VH.d(11:17) |  VL2.a(12:8)-CL(13:9){1}|CL(14:4){1}-VL.c(15:5)|CL(16:10){1}-VL.d(17:11)",
"scFV-IgG-scFV-scFV": "VL2.a(1:2)-L-VH2.a(2:1)-L-VH.a(3:12)-CH1(4:13){1}-H(5:18){2}-CH2(6:19)-CH3(7:20)-L-VH.b(8:9)-L-VL.b(9:8)-L-VH.c(10:11)-L-VL.c(11:10)|VL.a(12:3)-CL(13:4){1}|VL2.a(14:15)-L-VH2.a(15:14)-L-VH.a(16:25)-CH1(17:26){1}-H(18:5){2}-CH2(19:6)-CH3(20:7)-L-VH.b(21:22)-L-VL.b(22:21)-L-VH.c(23:24)-L-VL.c(24:23)|VL.a(25:16)-CL(26:17){1}",
"scFV-scFV-Fc":"VH.b(1:2)-L-VL.b(2:1)-L-VH.a(3:4)-L-VL.a(4:3)-CH2(5:7)-CH3(6:8)-L-CH2(7:5)-CH3(8:6)",
"Trimeric Fusion Protein":"X(8)[NOTE:FUSION]-X(9:14)[NOTE:FUSION]-CH1(10:15){1}-H(11:3){2}-CH2(12:4)-CH3(13:5)|X(14:9)[NOTE:FUSION]-CL(15:10){1}|VH.a(1:6)-CH1(2:7){1}-H(3:11){2}-CH2(4:12)-CH3(5:13)|VL.a(6:1)-CL(7:2){1}"
}

formats_keyslist= list(antibodyformats.keys())
for i in range(len(formats_keyslist)):
    Library.insert("end",formats_keyslist[i])
    if i %2==0:
        Library.itemconfig(i, bg='#D3D3D3')


Library.place(relx=0.01, rely = 0.425,relheight = 0.2, relwidth = 0.4)
Library.bind('<<ListboxSelect>>', items_selected)

#Library.place(frame1, relheight=0.05, relwidth=0.3)
###Results canvas

lower_frame = tk.Frame(root, bg = '#80c1ff', bd=5)
lower_frame.place(relx=0.45, rely=0.015, relwidth=0.55,relheight=0.93)
scrollbar = tk.Scrollbar(root)
scrollbar.pack( side = "right", fill = "y" )
lower_canvas = tk.Canvas(lower_frame, yscrollcommand = scrollbar.set)
lower_canvas.place(relheight=1,relwidth=1)
mm = MouseMover()
canvas_polygons = {}
canvas_labels   = {}
temp_label      = {}
custom_labels   = {}
# Bind mouse events to methods (could also be in the constructor)
lower_canvas.bind("<Button-1>", mm.select)
lower_canvas.bind("<B1-Motion>", mm.drag)
lower_canvas.bind("<ButtonRelease-1>", mm.release)
lower_canvas.bind("<Button-2>", mm.change_orientation)
startcoordinates = mm.select
newcoordinates = mm.drag

export_frame = tk.Frame(root, bg='#FF0000')
export_frame.place(relx=0.79, rely=0.945, relwidth=0.20,relheight=0.03)
template_file_button = tk.Button(export_frame, text = "Export template file", bg = "grey", font=40, command=lambda: Get_Template_File(lower_canvas))
template_file_button.place(relx=0, rely=0, relwidth=0.5,relheight=1)
Image_file_button = tk.Button(export_frame, text = "Export PNG", bg = "grey", font=40, command=lambda: save_as_png(lower_canvas))
Image_file_button.place(relx=0.5, rely=0, relwidth=0.5,relheight=1)







root.mainloop()
