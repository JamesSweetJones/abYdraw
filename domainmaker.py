#!/usr/bin/python
import re
import sys
import os
import tkinter as tk
from graphics import *

HEIGHT = 900
WIDTH  = 1200


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
    y       = re.sub("\s","",x)
    if "|" in y:
        splitx  = y.split("|")
    else:
        splitx = [y]
    if len(splitx) == 4:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a']
    elif len(splitx) == 2:
        chains = ['VH.b','VH.a']
    elif len(splitx) == 1:
        chains = ['VH.a']
    ADCs    = []
    VHa     = {}
    VHb     = {}
    VLa     = {}
    VLb     = {}
    fragment= {}
    Salt_bridges   = ""
    VHa_VLa_bonds  = ""
    VHb_VLb_bonds  = ""
    CH1a_CL1a_bonds=""
    CH1b_CL1b_bonds=""
    #get chains and watch out for ADCs
    #for i in splitx:
    #    if len(splitx) > 1:
    #        if i[0] != "X":
    #            chain = i.split("(")[0]
    #            chains.append(chain)
    #        elif i[0] == "X":
    #            chain = i.split("-")[1]
    #            chain = chain.split("(")[0]
    #            chains.append(chain)
    #    elif len(splitx) == 1:
    #        chains.append("fragment")

    for i in range(len(chains)):

        chain = splitx[i].split("-")
        if  "VH" in chain[0]:
            if chain[1] != "L":
                dict = str(re.sub("\.|\+|\_","",str(chains[i])))
            elif chain[1] == "L":
                if "VH" in chain[2]:
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                elif "VL" in chain[2]:
                    dict = str(re.sub("\.|\+|\_","",str(chains[i])))

        elif "VL" in chain[0]:
            if  chain[1] != "L":
                dict = str(re.sub("\.|\+|\_","",str(chains[i])))
            elif chain[1] == "L":
                if len(chain) > 2:
                    if "VH" in chain[2]:
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))
                    elif "VL" in chain[2]:
                        dict = str(re.sub("\.|\+|\_","",str(chains[i])))

        elif "X" in chain[0] and "VH" in chain[1]:
             dict = str(re.sub("\.|\+|\_","",str(chains[i])))
        #elif "VH" in chain[0] and
        for j in range(len(chain)):
            domain   =  re.findall("^CH[0-9][@+>_!*]\(.*\)\[.*\]|^CH[0-9]\(.*\)\[.*\]|^CH[0-9][@+>_!*]|^CL[0-9][@+>_!*]\(.*\)\[.*\]|^CL[0-9]\(.*\)\[.*\]|^CL[0-9][@+>_!*]|^CL[@+>_!*]|^CH[0-9]|^VL[1-9].[abcd]|^VL.[abcd]|^VH\.[abcd]|^VL[1-9]|^VH[1-9]\.[abcd]|^VH[1-9]|VH[+_*]\.[abcd]|VL[+_*]\.[abcd]|^CL[0-9]|^CL|^VH|^VL|^H|^X|^L", str(chain[j]))
            for i in range(len(domain)):
                if "*" in domain[i]:
                    domain = str(re.sub("\(.*\)\[MOD\:","",str(domain)))
                elif ("VH+") in domain[i] or ("VL+") in domain[i]:
                    domain = str(re.sub("\+","",domain[i]))
                    domain = domain+"+"
                elif ("VH_") in domain[i] or ("VL_") in domain[i]:
                    domain = str(re.sub("\_","",domain[i]))
                    domain = domain+"_"
            else:
                domain = str(re.sub("\[|\'|\]|\.","", str(domain)))


            #if chain[i] != fragment:
            #    if domain == ""
            ##Get Bond numbers
            if domain == "H" and Salt_bridges == "":
                Salt_bridges = re.findall("\{.*?\}", str(chain[j]))
                Salt_bridges = str(Salt_bridges)
                Salt_bridges = int(re.sub("\{|\'|\}|\[|\]","", Salt_bridges))

            elif domain == "VH":
                if dict == "VHa" and VHa_VLa_bonds == "":
                    VHa_VLa_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if VHa_VLa_bonds_match != []:
                        VHa_VLa_bonds = str(VHa_VLa_bonds_match)
                        VHa_VLa_bonds = int(re.sub("\{|\'|\}|\[|\]","", VHa_VLa_bonds))
                elif dict == "VHb" and VHb_VLb_bonds == "":
                    VHb_VLb_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if VHb_VLb_bonds_match != []:
                        VHb_VLb_bonds = str(VHb_VLb_bonds_match)
                        VHb_VLb_bonds = int(re.sub("\{|\'|\}|\[|\]","", VHb_VLb_bonds))
            elif domain == "CH1":
                if dict == "VHa" and CH1a_CL1a_bonds == "":
                    CH1a_CL1a_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if CH1a_CL1a_bonds_match != []:
                        CH1a_CL1a_bonds = str(CH1a_CL1a_bonds_match)
                        CH1a_CL1a_bonds = int(re.sub("\{|\'|\}|\[|\]","", CH1a_CL1a_bonds))
                elif dict == "VHb" and CH1b_CL1b_bonds == "":
                    CH1b_CL1b_bonds_match = re.findall("\{.*?\}", str(chain[j]))
                    if CH1b_CL1b_bonds_match != []:
                        CH1b_CL1b_bonds = str(CH1b_CL1b_bonds_match)
                        CH1b_CL1b_bonds = int(re.sub("\{|\'|\}|\[|\]","", CH1b_CL1b_bonds))

            location    = []
            locationstr =  re.findall("\((.*?)\)", str(chain[j]))
            locationstr =  str(re.sub("\[|\'|\]","", str(locationstr)))
            locationstr =  str(re.sub(":",",", str(locationstr)))
            locationstr = locationstr.split(",")

            if domain != "X" and domain !="L" and domain != "":
                for x in range(len(locationstr)):
                    location.append(int(locationstr[x]))
            elif domain == "X":
                location = re.findall("\[(.*?)\]", str(chain[j]))
            elif domain == "L":
                domain = "L"+str([j])
                location.append("")

            if dict == "VHa" and domain !="":
                VHa[domain] = location
            elif dict == "VHb" and domain !="":
                VHb[domain] = location
            elif dict == "VLa" and domain !="":
                VLa[domain] = location
            elif dict == "VLb" and domain !="":
                VLb[domain] = location
            elif dict == "fragment" and domain !="":
                fragment[domain] = location
            else:
                continue

    return(VHa,VLa,VHb,VLb,Salt_bridges,VHa_VLa_bonds,VHb_VLb_bonds,CH1a_CL1a_bonds,CH1b_CL1b_bonds,fragment)

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
    Salt_bridge_count     = chains_list[4]
    VHa_VLa_bond_count    = chains_list[5]
    VHb_VLb_bond_count    = chains_list[6]
    CH1a_CLa_bond_count   = chains_list[7]
    CH1b_CLb_bond_count   = chains_list[8]
    fragment              = chains_list[9]



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


    def domainmaker(startx,starty,righthanded,slant,VH,VL,X,mod,interaction):
        if VL == False and VH == False and mod == "" or X==True:
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
        elif VL == True:
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
        elif VH == True:
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
        elif VL == False and VH == False and mod == "@":
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
                sixthx     = fifthx + 40
            elif slant == False and righthanded == True:
                sixthx     = fifthx - 40
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx + 20
            elif righthanded == True:
                seventhx   =  firstx - 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]

        elif VL == False and VH == False and mod == ">" and interaction == "H_L":
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
                sixthx     = fifthx + 40
            elif slant == False and righthanded == True:
                sixthx     = fifthx - 40
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx - 20
            elif righthanded == True:
                seventhx   =  firstx + 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
        elif VL == False and VH == False and mod == ">" and interaction != "H_L":
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
                sixthx     = fifthx - 20
            elif slant == False and righthanded == True:
                sixthx     = fifthx + 20
            sixthy     = fifthy-37
            if righthanded == False:
                seventhx   =  firstx + 20
            elif righthanded == True:
                seventhx   =  firstx - 20
            seventhy   =  firsty
            coordinates = [firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy, sixthx, sixthy, seventhx, seventhy]
        elif VL == False and VH == False and X == True:
            firstx      = startx
            firsty      = starty

        top_bond    = [firstx,firsty]
        bottom_bond = [fourthx, fourthy]
        return(coordinates, bottom_bond,top_bond)



    def renderchains(dictionary,startx,starty):
        chain = []
        bonds = []
        slant = True
        righthanded = False
        Light_chain_check = False
        if dictionary == VHb_chain or dictionary == VLb_chain:
            Heavy_chain = VHb_chain
            Light_chain = VLb_chain
            righthanded = True
        elif dictionary == VHa_chain or dictionary == VLa_chain:
            Heavy_chain = VHa_chain
            Light_chain = VLa_chain
        if dictionary == VLa_chain or dictionary == VLb_chain:
            Light_chain_check = True


        for i in range(len(dictionary)):
            keyslist = list(dictionary.keys())
            V  = False
            VH = False
            VL = False
            X  = False
            print(keyslist[i])

            interaction = ""
            mod= ""
            if "VL" in keyslist[i]:
                VL = True
            elif "VH" in keyslist[i]:
                VH = True
            elif "X" in keyslist[i]:
                X  = True


            if "@" in keyslist[i]:
                mod = "@"
            elif ">" in keyslist[i]:
                mod = ">"
                if "H" in str(dictionary):
                    interaction = VH_VL_check(keyslist[i],dictionary,Light_chain_check,Heavy_chain,Light_chain,">","@")
            elif "+" in keyslist[i]:
                mod = "+"
            elif "_" in keyslist[i]:
                mod = "-"
            elif "*" in keyslist[i]:
                try:
                    mod = str(keyslist[i].split("*")[1])
                except IndexError:
                    mod = ""



            if i == 0:
                getcoordinates = domainmaker(startx,starty, righthanded,slant,VH,VL,X,mod,interaction)
                chain.append(getcoordinates[0])
                bottom_bond = getcoordinates[1]


            elif i > 0:
                previous_domain = keyslist[i-1]
                previous_chain  = chain[i-1]
                if keyslist[i-1] != "H" and "L[" not in keyslist[i-1] and "X" not in keyslist[i-1]:
                    previous_number = (dictionary.get(previous_domain)[0])+1
                elif keyslist[i-1] == "X" and i-1 == 0:
                    previous_number = (dictionary.get(keyslist[i])[0])
                    #previous_domain = [(dictionary.get(keyslist[i])[0])-1,0]
                    dictionary[keyslist[i-1]] = [previous_number,0]
                    previous_domain = (keyslist[i])
                if previous_domain == "H":
                    previous_domain = keyslist[i-2]
                    previous_chain  = chain[i-2]
                    previous_number = (dictionary.get(previous_domain)[0])+2
                if keyslist[i] != "X":
                    print(dictionary.get(keyslist[i])[0], previous_number, dictionary.get(keyslist[i])[0], previous_domain[1])
                if keyslist[i-1] == "CL":
                    slant = False
                if keyslist[i] == "H":
                    slant = False
                    chain.append([])
                    bonds.append([])
                elif "L[" in keyslist[i] :
                    chain.append([])
                    bonds.append([])
                elif keyslist[i] != "H" and "L[" not in keyslist[i-1] and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    if slant == True and righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,VH,VL,X,mod,interaction)
                    elif slant == True and righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,VH,VL,X,mod,interaction)
                    else:
                        getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,VH,VL,X,mod,interaction)

                    chain.append(getcoordinates[0])
                    top_bond = getcoordinates[2]
                    bonds.append(bottom_bond + top_bond)
                    bottom_bond = getcoordinates[1]





                elif "L[" in keyslist[i-1]:
                    previous_domain = keyslist[i-2]
                    if dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                        previous_chain = chain[i-2]
                        if slant == True and righthanded == True:
                            getcoordinates = domainmaker((previous_chain[6])-5,(previous_chain[7]+20),righthanded,slant,VH,VL,X,mod,interaction)
                        elif slant == True and righthanded == False:
                            getcoordinates = domainmaker((previous_chain[6])+5,(previous_chain[7]+20),righthanded,slant,VH,VL,X,mod,interaction)
                        else:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,VH,VL,X,mod,interaction)

                        chain.append(getcoordinates[0])
                        top_bond = getcoordinates[2]
                        bonds.append(bottom_bond + top_bond)
                        bottom_bond = getcoordinates[1]

                    elif dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] == (dictionary.get(previous_domain)[1]):
                        previous_chain = chain[i-2]
                        if "VH" in keyslist[i-2] and "VL" in keyslist[i]:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0]-50),(previous_chain[1]), righthanded,slant,VH,VL,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]+50),(previous_chain[1]), righthanded,slant,VH,VL,X,mod,interaction)
                        elif "VL" in keyslist[i-2] and "VH" in keyslist[i]:
                            if righthanded == False:
                                getcoordinates = domainmaker((previous_chain[0]+50),(previous_chain[1]), righthanded,slant,VH,VL,X,mod,interaction)
                            elif righthanded == True:
                                getcoordinates = domainmaker((previous_chain[0]-50),(previous_chain[1]), righthanded,slant,VH,VL,X,mod,interaction)
                        chain.append(getcoordinates[0])
                        top_bond = getcoordinates[2]
                        bonds.append(bottom_bond + top_bond)
                        bottom_bond = getcoordinates[1]
            #print(chain)







        chain = [x for x in chain if x != []]
        bonds = [x for x in bonds if x != []]
        return(chain,bonds)
    VHa_Domains_before_H = 0
    VHb_Domains_before_H = 0
    for i in range(len(VHa_chain)):
        keyslist = list(VHa_chain.keys())
        if keyslist[i] != "H" and "L[" not in keyslist[i]:
            VHa_Domains_before_H += 1
        elif keyslist[i] == "H":
            break
    for i in range(len(VHb_chain)):
        keyslist = list(VHb_chain.keys())
        if keyslist[i] != "H" and "L[" not in keyslist[i]:
            VHb_Domains_before_H += 1
        elif keyslist[i] == "H":
            break

    VLa_Domains_before_CL = 0
    VLb_Domains_before_CL = 0
    for i in range(len(VLa_chain)):
        keyslist = list(VLa_chain.keys())
        if  "CL" not in keyslist[i] and "L[" not in keyslist[i]:
            VLa_Domains_before_CL += 1
        elif "CL" in keyslist[i]:
            break
    for i in range(len(VLb_chain)):
        keyslist = list(VLb_chain.keys())
        if "CL" not in keyslist[i]  and "L[" not in keyslist[i]:
            VLb_Domains_before_CL += 1
        elif "CL" in keyslist[i]:
            break
    VLa_startx, VLa_starty = 75, 95
    VHa_startx, VHa_starty = 145, 95
    if VHa_Domains_before_H   == 2 and VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 95
        VHa_startx, VHa_starty = 145, 95
    elif VHa_Domains_before_H == 2 and VLa_Domains_before_CL == 2:
        VLa_startx, VLa_starty = 30,0
        VHa_startx, VHa_starty = 100,0
    elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 95
        VHa_startx, VHa_starty = 100,0
    elif VHa_Domains_before_H == 4 and VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty  = 75,95
        VHa_startx, VHa_starty  = 30,0
    if VLa_Domains_before_CL == 1:
        VLa_startx, VLa_starty = 75, 95

    VLb_startx, VLb_starty = 520, 95
    VHb_startx, VHb_starty = 450, 95
    if VHb_Domains_before_H   == 2 and VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 520, 95
        VHb_startx, VHb_starty = 450, 95
    elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 2:
        VLb_startx, VLb_starty = 560,0
        VHb_startx, VHb_starty = 480,0
    elif VHa_Domains_before_H == 3 and VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 525, 95
        VHb_startx, VHb_starty = 480,0
    elif VHa_Domains_before_H == 4 and VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty  = 520, 95
        VHb_startx, VHb_starty  = 560,0
    if VLa_Domains_before_CL == 1:
        VLb_startx, VLb_starty = 520, 95

    print(renderchains(VLa_chain, VLa_startx,VLa_starty)[0] + renderchains(VHa_chain, VHa_startx,VHa_starty)[0] + renderchains(VHb_chain, VHb_startx,VHb_starty)[0] + renderchains(VLb_chain, VLb_startx,VLb_starty)[0])
    print(renderchains(VLa_chain, VLa_startx,VLa_starty)[1] + renderchains(VHa_chain, VHa_startx,VHa_starty)[1] + renderchains(VHb_chain, VHb_startx,VHb_starty)[1] + renderchains(VLb_chain, VLb_startx,VLb_starty)[1])



entry = Get_input(sys.argv[1])
split_chains = Get_dictionaries(entry)
coordinates  = Check_interactions(split_chains)
