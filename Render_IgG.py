#!/usr/bin/python
import re
import sys
import os
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
    y       = re.sub("\s","",x)
    if "|" in y:
        splitx  = y.split("|")
    else:
        splitx = [y]
    if len(splitx) == 4:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a']
    elif len(splitx) == 2:
        chains  = ['VH.b','VH.a']
    elif len(splitx) == 1:
        chains  = ['VH.a']
    elif len(splitx) ==5:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1']
    elif len(splitx) ==6:
        chains  = ['VH.b', 'VL.b', 'VH.a', 'VL.a', 'fragment1', 'fragment2']


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
            if  "VH" in chain[0]:
                if chain[1] != "L":
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

            elif "X" in chain[0] and "VH" in chain[1]:
                 dict = str(re.sub("\.|\+|\_","",str(chains[i])))
        elif i >= 4:
            dict_number = i-3
            dict = str("fragment"+str(dict_number))


        for j in range(len(chain)):
            domain   =  re.findall("^CH[0-9][@+>_!*]\(.*\)\[.*\]|^CH[0-9]\(.*\)\[.*\]|^CH[0-9][@+>_!*]|^CL[0-9][@+>_!*]\(.*\)\[.*\]|^CL[0-9]\(.*\)\[.*\]|^CL[0-9][@+>_!*]|^CL[@+>_!*]|^CH[0-9]|^VL[1-9].[abcd]|^VL.[abcd]|^VH\.[abcd]|^VL[1-9]|^VH[1-9]\.[abcd]|^VH[1-9]|VH[@>+_*]\.[abcd]|VL[@>+_*]\.[abcd]|^CL[0-9]|^CL|^VH|^VL|^H\*\(.*\)\{.*\}\[.*\]|^H\(.*\)\{.*\}\[.*\]|^H|^X|^L", str(chain[j]))
            for i in range(len(domain)):
                mod = ""
                if "*" in domain[i]:
                    print(domain[i])
                    mod = str(re.sub("^.*\[|\].*","",str(domain[i])))
                    domain = str(re.sub("\*.*","",str(domain)))
                elif ("VH+") in domain[i] or ("VL+") in domain[i]:
                    domain = str(re.sub("\+","",domain[i]))
                    domain = domain+"+"
                elif ("VH_") in domain[i] or ("VL_") in domain[i]:
                    domain = str(re.sub("\_","",domain[i]))
                    domain = domain+"-"
            else:
                domain = str(re.sub("\[|\'|\]|\.","", str(domain)))


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

            if re.findall("\{.*?\}", str(chain[j])) != []:
                disulphide_bridges = re.findall("\{.*?\}", str(chain[j]))
                disulphide_bridges = re.sub("\{|\'|\}|\[|\]","", str(disulphide_bridges))
                disulphide_bridges = int(disulphide_bridges)
                location = [location,disulphide_bridges,mod]
            else:
                disulphide_bridges = 0
                location = [location, disulphide_bridges,mod]


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
                fragment4[domain] = location
            elif dict == "fragment4" and domain !="":
                fragment1[domain] = location
            else:
                continue

    return(VHa,VLa,VHb,VLb,chain_count,fragment1,fragment2,fragment3,fragment4)

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

    def innie_or_outie(chain,Light_chain_check, chain_count,Build_in,Build_out, equal_chain_lengths):
        innie_or_outie_list = []
        keyslist = list(chain.keys())

        for n in range(len(chain)):
            if "V" in keyslist[n]:
                default = ""
                if Light_chain_check == False:
                    default = "outie"
                elif Light_chain_check == True:
                    default = "innie"
                if n == 0:
                    try:
                        if "L[" in keyslist[n+1]:
                            if len(chain.get(keyslist[n])[0]) == 1:
                                innie_or_outie_list.append("Single_Fv_Chain")
                            else:
                                try:
                                    if  chain.get(keyslist[n])[0][0] == (chain.get(keyslist[n+2])[0][1]) and "L[" not in keyslist[n+3]:
                                        if Light_chain_check == False and Build_in == True:
                                            innie_or_outie_list.append("innie")
                                        elif Light_chain_check==False and Build_out== True:
                                            innie_or_outie_list.append("outie")
                                        elif Light_chain_check == True:
                                            innie_or_outie_list.append("innie")
                                    elif chain.get(keyslist[n])[0][0] == (chain.get(keyslist[n+2])[0][1]) and "L[" in keyslist[n+3]:
                                        if Light_chain_check == False and Build_out == True:
                                            innie_or_outie_list.append("outie")
                                        if Light_chain_check == False and Build_in == True:
                                            innie_or_outie_list.append("innie")

                                    elif chain.get(keyslist[n])[0][0] != (chain.get(keyslist[n+2])[0][1]):
                                        if chain_count == 2:
                                            if Light_chain_check == False:
                                                innie_or_outie_list.append("innie")
                                            elif Light_chain_check == True:
                                                innie_or_outie_list.append("innie")
                                        else:
                                            innie_or_outie_list.append(default)
                                except IndexError:
                                    innie_or_outie_list.append("innie")
                        elif "L[" not in keyslist[n+1]:
                            if len(chain.get(keyslist[n])[0]) == 1:
                                innie_or_outie_list.append("Single_Fv_Chain")
                            elif chain_count == 2 or chain_count == 1:
                                if Light_chain_check == False:
                                    innie_or_outie_list.append("innie")
                            else:
                                innie_or_outie_list.append(default)
                    except IndexError:
                        innie_or_outie_list.append(default)



                elif n > 0 and n+1 < len(chain):
                    if "L[" in keyslist[n-1]:
                        if len(chain.get(keyslist[n])[0]) == 1:
                            innie_or_outie_list.append("Single_Fv_Chain")
                        elif keyslist[n-2] == "X" and chain_count == 1:
                            innie_or_outie_list.append("innie")
                        elif chain.get(keyslist[n-2])[0][0] == (chain.get(keyslist[n])[0][1]):

                            if innie_or_outie_list[-2] == "outie":
                                innie_or_outie_list.append("innie")
                            elif innie_or_outie_list[-2] == "innie":
                                innie_or_outie_list.append("outie")
                            else:
                                innie_or_outie_list.append(default)

                        elif  chain.get(keyslist[n-2])[0][0] != (chain.get(keyslist[n])[0][1]):
                            if "CL" in keyslist[n-2] and Light_chain_check == True and Build_out == True:
                                innie_or_outie_list.append("outie")
                            elif innie_or_outie_list[-2] == "outie":
                                innie_or_outie_list.append("outie")
                            elif innie_or_outie_list[-2] == "innie" and (chain_count==1 or chain_count==2):
                                innie_or_outie_list.append("innie")
                            elif innie_or_outie_list[-2] == "innie" and chain_count==4 and Light_chain_check == False:
                                if Build_in==True:
                                    innie_or_outie_list.append("innie")
                                elif Build_out==True:
                                    innie_or_outie_list.append("outie")

                            else:
                                innie_or_outie_list.append(default)
                        else:
                            innie_or_outie_list.append(default)

                    else:
                        innie_or_outie_list.append(default)

                elif n+1 == len(chain):
                    if len(chain.get(keyslist[n])[0]) == 1:
                        innie_or_outie_list.append("Single_Fv_Chain")
                    elif chain.get(keyslist[n-2])[0][0] == (chain.get(keyslist[n])[0][1]):
                        if innie_or_outie_list[-2] == "outie":
                            innie_or_outie_list.append("innie")
                        elif innie_or_outie_list[-2] == "innie":
                            innie_or_outie_list.append("outie")
                        else:
                            innie_or_outie_list.append(default)
                    elif "X" in str(keyslist):
                        innie_or_outie_list.append("outie")

                    else:
                        innie_or_outie_list.append("innie")

            elif "V" not in keyslist[n]:
                innie_or_outie_list.append("constant")


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
                sixthx     = fifthx + 20
            elif slant == False and righthanded == True:
                sixthx     = fifthx - 20
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

            coordinates=[firstx , firsty , secondx , secondy, thirdx, thirdy, fourthx, fourthy, fifthx, fifthy,sixthx,sixthy,seventhx,seventhy,eighthx,eighthy,ninthx,ninthy]
        elif mod == "H":
            firstx, firsty = startx,starty
            secondx, secondy=startx,starty
            thirdx, thirdy=startx,starty
            fourthx, fourthy=startx,starty
            fifthx, fifthy=startx,starty

            coordinates = []
        if slant == True or mod == "Leucine" or previous_H == True :
            top_bond    = [firstx,firsty]
        else:
            top_bond    = [firstx,firsty+20]

        if mod=="Leucine":
            bottom_bond = [ninthx,ninthy]
        else:
            bottom_bond = [fourthx, fourthy]

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



    def renderchains(dictionary,startx,starty):
        chain                   = []
        coordinates_list_heavy_a= []
        coordinates_list_light_a= []
        coordinates_list_heavy_b= []
        coordinates_list_light_b= []
        coordinates_list_heavy_c= []
        coordinates_list_light_c= []
        coordinates_list_heavy_d= []
        coordinates_list_light_d= []
        bonds = []
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
        ADCs            =[]
        keyslist = list(dictionary.keys())


        if chain_count >= 4:
            slant = True
        elif chain_count==2 and "H" in dictionary and len(keyslist) > 4:
            slant = True
        else:
            slant = False
        before_H = True

        if chain_count == 4:
            if dictionary == VHa_chain or dictionary == VHb_chain:
                keyslist = list(dictionary.keys())
                try:
                    if "L[" in keyslist[1] and dictionary.get(keyslist[0])[0][0] == (dictionary.get(keyslist[2])[0][1]):
                        Build_out = True
                        Build_in  = False
                    else:
                        Build_in = True
                        Build_out = False
                except IndexError:
                    Build_in = True
                    Build_out = False
            else:
                Build_in = True
                Build_out = False

        else:
            Build_in = True
            Build_out = False
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
        innie_or_outie_list = innie_or_outie(dictionary, Light_chain_check, chain_count,Build_in,Build_out,equal_chain_lengths)
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
            Build_up=False
            Build_down=True


            if dictionary.get(keyslist[i])[2] != "":
                print("HELLO THERE")
                note_label = str(dictionary.get(keyslist[i])[2])
                Notes.append(str(keyslist[i]+" "+note_label))
                if len(Notes_positions) == 0:
                    Notes_positions.append([200,600])
                elif len(Notes_positions) > 0:
                    XY = (Notes_positions[-1][1])+20
                    Notes_positions.append([200,XY])

            print(mod_label, mod)



            try:
                location = dictionary.get(keyslist[i])[0][0]
                disulphide_bridge_count = int(dictionary.get(keyslist[i])[1])
                location    = dictionary.get(keyslist[i])[0]
                dictionary[keyslist[i]] = location
                if keyslist[i] == "H":
                    H_disulphide_bridge_count = disulphide_bridge_count

            except:
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




            elif "*" in keyslist[i]:
                try:
                    mod_label = str(keyslist[i].split("*")[1])
                    Notes.append(mod_label)
                    if len(Notes_positions) == 0:
                        Notes_positions.append([200,600])
                    elif len(Notes_positions) > 0:
                        XY = (Notes_positions[-1][1])+20
                        Notes_positions.append([200,XY])
                except IndexError:
                    mod = ""


            Domain_name = str(re.sub("\@|\>|\<","",str(keyslist[i])))
            #print(keyslist[i], i, len(dictionary))
            print(keyslist[i])


            if i == 0:
                print("checkpoint1")
                getcoordinates = domainmaker(startx,starty, righthanded,slant,V,direction,X,mod,interaction,previous_H)

            elif i > 0:
                previous_domain = keyslist[i-1]
                previous_chain  = chain[i-1]
                if keyslist[i-1] != "H" and "L[" not in keyslist[i-1] and "X" not in keyslist[i-1]:
                    previous_number = (dictionary.get(previous_domain)[0])+1
                elif "X" in keyslist[i-1]  and i-1 == 0:
                    previous_number = (dictionary.get(keyslist[i])[0])
                    dictionary[keyslist[i-1]] = [previous_number,0]
                    previous_domain = (keyslist[i])
                elif "X" in keyslist[i-1] and chain_count > 2:
                    slant=False


                if previous_domain == "H":
                    previous_domain = keyslist[i-2]
                    previous_chain  = chain[i-2]
                    print(dictionary.get(previous_domain))
                    previous_number = int(dictionary.get(previous_domain)[0])+2


                if dictionary == VLa_chain or dictionary == VLb_chain:
                    if keyslist[i-1] == "CL":
                        slant = False

                if keyslist[i] == "H":
                    before_H = False


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



                elif keyslist[i-1] == "H" and keyslist[i]== "X":
                    print("checkpoint4")
                    if righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)



                elif "X" in keyslist[i] and keyslist[i-1] != "H":
                    print("checkpoint5")
                    if "L[" in keyslist[i-1]:
                        previous_chain = chain[i-2]
                    if "C" in keyslist[i-1]:
                        if Build_in == True:
                            Build_in = False
                            Build_out = True
                    if mod !="Leucine":

                        if righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker((previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker((previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker((previous_chain[0]-98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker((previous_chain[0]+98),(previous_chain[1]),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        Build_in = True
                        Build_out = False
                    elif mod =="Leucine":
                        if righthanded == False :
                            getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True :
                            getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif slant == False:
                            getcoordinates = domainmaker((previous_chain[6]),(previous_chain[7]+20),righthanded,slant,V,direction,X,mod,interaction,previous_H)



                elif keyslist[i-1] == "H" and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    print("checkpoint6")
                    previous_H = True
                    if righthanded == True:
                        getcoordinates = domainmaker((previous_chain[6]-20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    elif righthanded == False:
                        getcoordinates = domainmaker((previous_chain[6]+20),(previous_chain[7]+40),righthanded,slant,V,direction,X,mod,interaction,previous_H)


                elif "L[" not in keyslist[i-1] and "X" and len(dictionary.get(keyslist[i-1])) ==1:
                    print("checkpoint7")
                    if "L[" in keyslist[i]:
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

                elif keyslist[i] != "H" and "L[" not in keyslist[i-1] and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]) and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[1] == (dictionary.get(keyslist[i-1])[1]+2):
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


                elif keyslist[i] != "H" and "L[" not in keyslist[i-1] and dictionary.get(keyslist[i])[0] == previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
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
                elif "L[" in keyslist[i-1]:
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
                    elif len(dictionary.get(keyslist[i-2])) ==1 and keyslist[i-2] != "X":
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

                        if  righthanded == False and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == False and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif  righthanded == True and Build_in == True and Build_out == False:
                            getcoordinates = domainmaker(to_joinx-60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif righthanded == True and Build_in == False and Build_out == True:
                            getcoordinates = domainmaker(to_joinx+60,to_joiny, righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        all_list = list(All_positions_and_chains.keys())
                        get = (All_positions_and_chains.get(all_list[-1]))
                        yprev = get[0][1]

                        if dictionary.get(keyslist[i-2])[0] != dictionary.get(keyslist[i])[1] and to_joiny < yprev:
                            print("checkpoint13")
                            Build_up=True
                            Build_down=False

##ADCs
                    elif keyslist[i-2] == "X":
                        print("checkpoint14")
                        if chain_count ==1:
                            if Build_in == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker((previous_chain[0]+100),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif righthanded == True:
                                    getcoordinates = domainmaker((previous_chain[0]-100),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif Build_out == True:
                                if righthanded == False:
                                    getcoordinates = domainmaker((previous_chain[0]-100),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                                elif righthanded == True:
                                    getcoordinates = domainmaker((previous_chain[0]+100),(previous_chain[1]), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                        elif chain_count > 1:
                            if righthanded == False and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker((previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == False and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker((previous_chain[6]-50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True and Build_in == True and Build_out == False:
                                getcoordinates = domainmaker((previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
                            elif righthanded == True and Build_in == False and Build_out == True:
                                getcoordinates = domainmaker((previous_chain[6]+50),(previous_chain[7]+30),righthanded,slant,V,direction,X,mod,interaction,previous_H)
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
                                if "L[" in keyslist[i+1]:
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


                        if (in_out_counter%2)!=0:
                            Build_in  = False
                            Build_out = True
                        elif (in_out_counter%2)==0:
                            Build_in  = True
                            Build_out = False
##Build diagonally
                    #elif dictionary.get(keyslist[i])[0] != previous_number and dictionary.get(keyslist[i])[0] != (dictionary.get(previous_domain)[1]):
                    #    if righthanded  == True:
                    #        righthanded = False
                    #    elif righthanded == False:
                    #        righthanded = True
                    #    if Build_in == True:
                    #        if righthanded == False:
                    #            getcoordinates = domainmaker((previous_chain[6]-100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    #        elif righthanded == True:
                    #            getcoordinates = domainmaker((previous_chain[6]+100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    #    if Build_in == False:
                    #        if righthanded == False:
                    #            getcoordinates = domainmaker((previous_chain[6]+100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction,previous_H)
                    #        elif righthanded == True:
                    #            getcoordinates = domainmaker((previous_chain[6]-100),(previous_chain[7]+20), righthanded,slant,V,direction,X,mod,interaction,previous_H)

##H disulphides
            if keyslist[i] == "H" or keyslist[i] == "H*"  :
                if H_disulphide_bridge_count > 0 and Build_up == False:
                    bottomx=bottom_bond[0]
                    bottomy=bottom_bond[1]
                    if slant == True and righthanded == False:
                        topx = bottomx+20
                    elif slant == True and righthanded == True:
                        topx = bottomx-20
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






                slant=False

##append coordinates to chains
            if keyslist[i] == "H" or "L[" in keyslist[i]:
                chain.append([])
                bonds.append([])
            else:
                chain.append(getcoordinates[0])

            if i > 0 and Build_down==True:
                print("checkpoint19")
                top_bond = getcoordinates[2]
                bonds.append(bottom_bond + top_bond)
                if mod=="Leucine":
                    bonds.append(getcoordinates[0])
            elif i > 0 and Build_up==True :
                print("checkpoint20")
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
                        if H_disulphide_bridge_count > 10:
                            if noted_Hbonds==False:
                                Notes.append((H_disulphide_bridge_count, "hinge disulphide bonds"))
                                noted_Hbonds==True
                                if len(Notes_positions) == 0:
                                    Notes_positions.append([200,600])
                                elif len(Notes_positions) > 0:
                                    XY = (Notes_positions[-1][1])+20
                                    Notes_positions.append([200,XY])



            bottom_bond = getcoordinates[1]
            Build_up_downlist.append(Build_up)





###append domains to right list
            if keyslist[i] != "H" and "L[" not in keyslist[i] and "X" not in keyslist[i]:
                if i == 0:
                    if "a" in keyslist[i]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:
                            coordinates_list_light_a.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_a.append(getcoordinates[0])
                    elif "b" in keyslist[i]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:

                            coordinates_list_light_b.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_b.append(getcoordinates[0])
                    elif "c" in keyslist[i]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:
                            coordinates_list_light_c.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_c.append(getcoordinates[0])
                    elif "d" in keyslist[i]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:
                            coordinates_list_light_d.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_d.append(getcoordinates[0])
                    elif "CL" in keyslist[i]:
                        if dictionary == VHa_chain or dictionary == VLa_chain:
                            coordinates_list_light_a.append(getcoordinates[0])
                        elif dictionary == VHb_chain or dictionary == VLb_chain:
                            coordinates_list_light_b.append(getcoordinates[0])
                        elif dictionary == fragment1:
                            if "a"  in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_a.append(getcoordinates[0])
                            elif "a"  in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_a.append(getcoordinates[0])
                            elif "b"  in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_b.append(getcoordinates[0])
                            elif "b"  in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_b.append(getcoordinates[0])
                            elif "c" in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_c.append(getcoordinates[0])
                            elif "c" in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_c.append(getcoordinates[0])
                            elif "d" in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_d.append(getcoordinates[0])
                            elif "d" in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_light_d.append(getcoordinates[0])
                        elif dictionary == fragment2:
                            if "a"  in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_a.append(getcoordinates[0])
                            elif "a"  in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_a.append(getcoordinates[0])
                            elif "b"  in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_b.append(getcoordinates[0])
                            elif "b"  in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_b.append(getcoordinates[0])
                            elif "c" in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_c.append(getcoordinates[0])
                            elif "c" in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_c.append(getcoordinates[0])
                            elif "d" in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_d.append(getcoordinates[0])
                            elif "d" in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_d.append(getcoordinates[0])
                    elif ("a" and "b" and "c" and "d") not in keyslist[i]:

                        if "VH" in keyslist[i]:
                            coordinates_list_heavy_a.append(getcoordinates[0])
                        elif "VL" in keyslist[i]:
                            coordinates_list_light_a.append(getcoordinates[0])
                        elif "X" in keyslist[i] and mod != "Leucine":
                            if "a" in str(keyslist):
                                coordinates_list_heavy_a.append(getcoordinates[0])
                            elif "b" in str(keyslist):
                                coordinates_list_heavy_b.append(getcoordinates[0])
                        else:
                            if "a"  in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_a.append(getcoordinates[0])
                            elif "a"  in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_a.append(getcoordinates[0])
                            elif "b"  in str(keyslist) and "H" in keyslist[i]:
                                coordinates_list_heavy_b.append(getcoordinates[0])
                            elif "b"  in str(keyslist) and "L" in keyslist[i]:
                                coordinates_list_light_b.append(getcoordinates[0])
                            else:
                                coordinates_list_heavy_a.append(getcoordinates[0])

                    else:
                        if ("a")  in str(keyslist) and "H" in str(keyslist[i]):
                            coordinates_list_heavy_a.append(getcoordinates[0])
                        elif ("a")  in str(keyslist) and "L" in str(keyslist[i]):
                            coordinates_list_light_a.append(getcoordinates[0])
                        elif ("b")  in str(keyslist) and "H" in str(keyslist[i]):
                            coordinates_list_heavy_b.append(getcoordinates[0])
                        elif ("b")  in str(keyslist) and "L" in str(keyslist[i]):
                            coordinates_list_light_b.append(getcoordinates[0])

                elif i > 0 and mod != "Leucine":
                    if "a" in keyslist[i] or "a" in keyslist[i-1]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:
                            coordinates_list_light_a.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_a.append(getcoordinates[0])
                    elif "b" in keyslist[i] or "b" in keyslist[i-1]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:

                            coordinates_list_light_b.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_b.append(getcoordinates[0])
                    elif "c" in keyslist[i]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:
                            coordinates_list_light_c.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_c.append(getcoordinates[0])
                    elif "d" in keyslist[i]:
                        if "VL" in keyslist[i] or "CL" in keyslist[i]:
                            coordinates_list_light_d.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_d.append(getcoordinates[0])
                    elif "CL" in keyslist[i]:
                        if "a" in str(keyslist):
                            coordinates_list_light_a.append(getcoordinates[0])
                        elif "b" in str(keyslist):
                            coordinates_list_light_b.append(getcoordinates[0])
                        else:
                            coordinates_list_light_a.append(getcoordinates[0])
                    elif  dictionary == VHa_chain or dictionary == VLa_chain:

                        coordinates_list_heavy_a.append(getcoordinates[0])
                    elif dictionary == VHb_chain or dictionary == VLb_chain:
                        if "b" in str(keyslist):
                            coordinates_list_heavy_b.append(getcoordinates[0])
                        else:
                            coordinates_list_heavy_a.append(getcoordinates[0])
                    elif dictionary == fragment1 or dictionary == fragment2:
                        if "a" in keyslist[i-1] or "a" in keyslist[i-2] :
                            coordinates_list_heavy_a.append(getcoordinates[0])
                        elif "b" in keyslist[i-1]or "b" in keyslist[i-2] :
                            coordinates_list_heavy_b.append(getcoordinates[0])
                        elif "c" in keyslist[i-1] or "c" in keyslist[i-2] :
                            coordinates_list_heavy_c.append(getcoordinates[0])
                        elif "d" in keyslist[i-1]or "d" in keyslist[i-2] :
                            coordinates_list_heavy_d.append(getcoordinates[0])





##Get labels and positions
            if keyslist[i] != "H" and "L[" not in keyslist[i] and "X" not in keyslist[i]:
                Label_Locations = getcoordinates[3]
                text_coordinates.append([Label_Locations])
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [getcoordinates[0], righthanded]
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
                text = dictionary.get(keyslist[i])
                Location_Text.append(str(text)+mod_label)
                Domain_Text.append(str(Domain_name)+mod_label)
                All_positions_and_chains[Domain_name] = [getcoordinates[0], righthanded]
                Notes.append((Domain_name, text))
                if len(Notes_positions) == 0:
                    Notes_positions.append([200,600])
                elif len(Notes_positions) > 0:
                    XY = (Notes_positions[-1][1])+20
                    Notes_positions.append([200,XY])

##Sort extra disulphide bridges
            if disulphide_bridge_count > 0 and keyslist[i] !="H" and keyslist[i] != "X":
                bottomx=bottom_bond[0]
                bottomy=bottom_bond[1]
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

            elif keyslist[i]=="H":
                print("OK OWWW")
                if i+1 != len(dictionary):
                    if Extra_bond == True:
                        Label_Locations = [getcoordinates[3][0],getcoordinates[3][1]+80]
                    else:
                        Label_Locations = [getcoordinates[3][0],getcoordinates[3][1]+60]
                elif i+1 == len(dictionary):
                    if Extra_bond == True and righthanded == False:
                        Label_Locations = [getcoordinates[3][0]-60,getcoordinates[3][1]]
                    elif Extra_bond == True and righthanded == True:
                        Label_Locations = [getcoordinates[3][0]+60,getcoordinates[3][1]]
                    elif righthanded== False:
                        Label_Locations = [getcoordinates[3][0]-60,getcoordinates[3][1]-17]
                    elif righthanded== True:
                        Label_Locations = [getcoordinates[3][0]+60,getcoordinates[3][1]-17]
                text_coordinates.append(Label_Locations)
                text = dictionary.get(keyslist[i])[0]
                Location_Text.append(str(text)+mod_label)
                All_positions_and_chains[text] = [Label_Locations, righthanded]
                Domain_Text.append(str(Domain_name)+mod_label)




        chain    = [x for x in chain if x != []]
        allbonds = bonds
        allbonds = [x for x in bonds if x != []]
        #allbonds = bonds
        return(coordinates_list_heavy_a,coordinates_list_light_a,coordinates_list_heavy_b,coordinates_list_light_b,coordinates_list_heavy_c,coordinates_list_light_c,coordinates_list_heavy_d,coordinates_list_light_d,allbonds,Location_Text,text_coordinates, disulphidebridge1, disulphidebridge2,disulphidebridge3,disulphidebridge4,disulphidebridge5,completed_disulphidebridges,Domain_Text,Notes,Notes_positions, arcs_left,arcs_right, arcs_left_slant, arcs_right_slant,ADCs)

    VLa_Domains_before_CL = 0
    CLa_count             = 0
    VLb_Domains_before_CL = 0
    CLb_count             = 0
    for i in range(len(VLa_chain)):
        keyslistlight_a = list(VLa_chain.keys())
        if  "CL" not in keyslistlight_a[i] and "L[" not in keyslistlight_a[i]:
            VLa_Domains_before_CL += 1
        elif "CL" in keyslistlight_a[i]:
            CLa_count +=1
            break
    for i in range(len(VLb_chain)):
        keyslistlight_b = list(VLb_chain.keys())
        if "CL" not in keyslistlight_b[i]  and "L[" not in keyslistlight_b[i]:
            VLb_Domains_before_CL += 1
        elif "CL" in keyslistlight_b[i]:
            CLb_count +=1
            break
    if CLa_count == 0 and VLa_Domains_before_CL > 1:
        VLa_startx, VLa_starty = 155, 110
    elif VLa_Domains_before_CL == 1 and len(VLa_chain)==1:
        VLa_startx, VLa_starty = 200,205
    elif VLa_Domains_before_CL == 2 and "L[" in keyslistlight_a[1] and (VLa_chain.get(keyslistlight_a[0])[0][0]) == (VLa_chain.get(keyslistlight_a[2])[0][1]):
        VLa_startx, VLa_starty = 95,110
    elif VLa_Domains_before_CL == 2 and  "L[" in keyslistlight_a[1] and (VLa_chain.get(keyslistlight_a[0])[0][0]) != (VLa_chain.get(keyslistlight_a[2])[0][1]):
        VLa_startx, VLa_starty = 105, 15

    else:
        VLa_startx, VLa_starty = 155, 110

    if CLb_count == 0 and VLb_Domains_before_CL > 1:
        VLb_startx, VLb_starty = 560, 110
    elif VLb_Domains_before_CL == 1 and len(VLb_chain)==1:
        VLb_startx, VLb_starty = 515,205
    elif VLb_Domains_before_CL == 2 and "L[" in keyslistlight_b[1] and (VLb_chain.get(keyslistlight_b[0])[0][0]) == (VLb_chain.get(keyslistlight_b[2])[0][1]):
        VLb_startx, VLb_starty = 610,110
    elif VLb_Domains_before_CL == 2 and "L[" in keyslistlight_b[1] and (VLb_chain.get(keyslistlight_b[0])[0][0]) != (VLb_chain.get(keyslistlight_b[2])[0][1]):
        VLb_startx, VLb_starty = 620, 15
    else:
        VLb_startx, VLb_starty = 560, 110



    VHa_Domains_before_H = 0
    VHb_Domains_before_H = 0
    keyslista = list(VHa_chain.keys())
    for i in range(len(VHa_chain)):
        if keyslista[i] != "H" and "L[" not in keyslista[i]:
            VHa_Domains_before_H += 1
        elif keyslista[i] == "H":
            break
    keyslistb = list(VHb_chain.keys())
    for i in range(len(VHb_chain)):
        if keyslistb[i] != "H" and "L[" not in keyslistb[i]:
            VHb_Domains_before_H += 1
        elif keyslistb[i] == "H":
            break
    if chain_count == 2:
        if "L[" in keyslista[1] and  VHa_chain.get(keyslista[0])[0][0] == (VHa_chain.get(keyslista[2])[0][1]) and "H" not in keyslista:
            VHa_startx, VHa_starty = 150,100
        elif "H" in keyslista:
            interactor = VHa_chain.get(keyslista[0])[0][1]
            interactor_found= False
            for i in range(1,len(keyslista)):
                try:
                    interaction = VHa_chain.get(keyslista[i])[0][0]
                    if interactor == interaction:
                        interactor_found=True
                        if VHa_Domains_before_H  < 3:
                            VHa_startx, VHa_starty = 200,205

                        else:
                            VHa_startx, VHa_starty = 150, 110
                except IndexError:
                    continue
            if interactor_found==False:
                VHa_startx, VHa_starty = 155,100
        else:
            VHa_startx, VHa_starty = 250, 100



    elif chain_count ==1:
        VHa_startx, VHa_starty = 350,100
        VHb_startx, VHb_starty = 500, 110
    elif chain_count >=4:
        if VHa_Domains_before_H == 1:
            VHa_startx, VHa_starty = 260, 205
        elif "X" in str(keyslista) and "X" in keyslista[0]:
            VHa_startx, VHa_starty = 215, 15
        elif "X" in str(keyslista) and "X" not in keyslista[0]:
            VHa_startx, VHa_starty = 215, 110
        elif VHa_Domains_before_H == 2:
            VHa_startx, VHa_starty = 215, 110
        elif VHa_Domains_before_H == 3 and "L[" in keyslista[1] and VHa_chain.get(keyslista[0])[0][0] == (VHa_chain.get(keyslista[2])[0][1]):
            VHa_startx, VHa_starty = 275, 110
        elif VHa_Domains_before_H == 3 and "L[" in keyslista[1]:
            VHa_startx, VHa_starty = 165, 15

        elif VHa_Domains_before_H == 3 and "L[" not in keyslista[1]:
            VHa_startx, VHa_starty = 215, 15
        elif VHa_Domains_before_H == 4 and "L[" in keyslista[1]:
            VHa_startx, VHa_starty = 225, 15
        elif VHa_Domains_before_H == 4 and "L[" not in keyslista[1]:
            VHa_startx, VHa_starty = 165, 15


    if chain_count == 2:
        if "L[" in keyslistb[1] and  VHb_chain.get(keyslistb[0])[0][0] == (VHb_chain.get(keyslistb[2])[0][1]) and "H" not in keyslistb:
            VHb_startx, VHb_starty = 400,100
        elif "H" in keyslistb:
            interactor = VHb_chain.get(keyslistb[0])[0][1]
            interactor_found= False
            for i in range(1,len(keyslistb)):
                try:
                    interaction = VHb_chain.get(keyslistb[i])[0][0]
                    if interactor == interaction:
                        interactor_found=True
                        if VHb_Domains_before_H  < 3:
                            VHb_startx, VHb_starty = 515,205

                        else:
                            VHb_startx, VHb_starty = 560, 110
                except IndexError:
                    continue
            if interactor_found==False:
                VHb_startx, VHb_starty = 400,100
        else:
            VHb_startx, VHb_starty = 350,100
    elif chain_count >=4:
        if VHb_Domains_before_H ==1:
            VHb_startx, VHb_starty = 455, 205
        elif "X" in str(keyslistb) and "X" not in keyslistb[0]:
            VHb_startx, VHb_starty = 500, 110
        elif "X" in str(keyslistb) and "X" in keyslistb[0]:
            VHb_startx, VHb_starty = 500, 15
        elif VHb_Domains_before_H ==2:
            VHb_startx, VHb_starty = 500, 110
        elif VHb_Domains_before_H == 3 and "L[" in keyslistb[1] and VHb_chain.get(keyslistb[0])[0][0] == (VHb_chain.get(keyslistb[2])[0][1]):
            VHb_startx, VHb_starty = 430, 110
        elif VHb_Domains_before_H ==3 and "L[" in keyslistb[1]:
            VHb_startx, VHb_starty = 550, 15
        elif VHb_Domains_before_H ==3 and not "L[" in keyslistb[1]:
            VHb_startx, VHb_starty = 500, 15
        elif VHb_Domains_before_H == 4 and "L[" in keyslistb[1]:
            VHb_startx, VHb_starty = 490, 15
        elif VHb_Domains_before_H == 4 and "L[" not in keyslistb[1]:
            VHb_startx, VHb_starty = 550, 15


    if chain_count ==2:
        VHb_startx, VHb_starty = VHb_startx, VHb_starty
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VLb_stats = renderchains(VLb_chain,VLb_startx, VLb_starty)
        VHa_list = list(VHa_chain.keys())

        try:
            VHa_inter = VHa_chain.get(VHa_list[0])[0][1]
            if len(VHa_chain) == len(VHb_chain) and VHa_chain.get(keyslista[0])[0][0] != VHa_chain.get(keyslista[2])[0][1]:
                VHa_start = find_the_fragment(VHa_inter,All_positions_and_chains)
                if VHa_start is not None:
                    righthanded = VHa_start[1]
                    if righthanded == True:
                        VHa_startx=VHa_start[0][0]-50
                        VHa_starty=VHa_start[0][1]
                    elif righthanded==False:
                        VHa_startx=VHa_start[0][0]-50
                        VHa_starty=VHa_start[0][1]
                else:
                    VHa_startx, VHa_starty= VHa_startx, VHa_starty
            else:
                VHa_startx, VHa_starty= VHa_startx, VHa_starty
        except IndexError:
            VHa_startx, VHa_starty= VHa_startx, VHa_starty

        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)
        VLa_stats = renderchains(VLa_chain,VLa_startx, VLa_starty)




    elif chain_count ==1 or chain_count >= 4:
        VHa_stats = renderchains(VHa_chain,VHa_startx, VHa_starty)
        VLa_stats = renderchains(VLa_chain,VLa_startx, VLa_starty)
        VHb_stats = renderchains(VHb_chain,VHb_startx, VHb_starty)
        VLb_stats = renderchains(VLb_chain,VLb_startx, VLb_starty)
    frag1_startx,frag1_starty,frag2_startx,frag2_starty,frag3_startx,frag3_starty,frag4_startx,frag4_starty = 0,0,0,0,0,0,0,0

    if fragment1 != {}:
        fragment1_list = list(fragment1.keys())
        fragment_inter = fragment1.get(fragment1_list[0])[0][1]
        frag1_start = find_the_fragment(fragment_inter,All_positions_and_chains)
        righthanded = frag1_start[1]
        if righthanded == True:
            frag1_startx=frag1_start[0][0]+50
            frag1_starty=frag1_start[0][1]
        elif righthanded==False:
            frag1_startx=frag1_start[0][0]-50
            frag1_starty=frag1_start[0][1]
    if fragment2 != {}:
        fragment2_list = list(fragment2.keys())
        fragment_inter = fragment2.get(fragment2_list[0])[0][1]
        frag2_start = find_the_fragment(fragment_inter,All_positions_and_chains)
        righthanded = frag2_start[1]
        if righthanded == True:
            frag2_startx=frag2_start[0][0]+50
            frag2_starty=frag2_start[0][1]
        elif righthanded==False:
            frag2_startx=frag2_start[0][0]-50
            frag2_starty=frag2_start[0][1]
    frag1_stat= renderchains(fragment1,frag1_startx,frag1_starty)
    frag2_stat= renderchains(fragment2,frag2_startx,frag2_starty)
    frag3_stat= renderchains(fragment3,frag3_startx,frag3_starty)
    frag4_stat= renderchains(fragment4,frag4_startx,frag4_starty)





    Heavy_Domains_a     = VHa_stats[0] + VLa_stats[0] + VHb_stats[0] + VLb_stats[0] + frag1_stat[0] + frag2_stat[0] + frag3_stat[0] + frag4_stat[0]
    Light_Domains_a     = VHa_stats[1] + VLa_stats[1] + VHb_stats[1] + VLb_stats[1] + frag1_stat[1] + frag2_stat[1] + frag3_stat[1] + frag4_stat[1]
    Heavy_Domains_b     = VHa_stats[2] + VLa_stats[2] + VHb_stats[2] + VLb_stats[2] + frag1_stat[2] + frag2_stat[2] + frag3_stat[2] + frag4_stat[2]
    Light_Domains_b     = VHa_stats[3] + VLa_stats[3] + VHb_stats[3] + VLb_stats[3] + frag1_stat[3] + frag2_stat[3] + frag3_stat[3] + frag4_stat[3]
    Heavy_Domains_c     = VHa_stats[4] + VLa_stats[4] + VHb_stats[4] + VLb_stats[4] + frag1_stat[4] + frag2_stat[4] + frag3_stat[4] + frag4_stat[4]
    Light_Domains_c     = VHa_stats[5] + VLa_stats[5] + VHb_stats[5] + VLb_stats[5] + frag1_stat[5] + frag2_stat[5] + frag3_stat[5] + frag4_stat[5]
    Heavy_Domains_d     = VHa_stats[6] + VLa_stats[6] + VHb_stats[6] + VLb_stats[6] + frag1_stat[6] + frag2_stat[6] + frag3_stat[6] + frag4_stat[6]
    Light_Domains_d     = VHa_stats[7] + VLa_stats[7] + VHb_stats[7] + VLb_stats[7] + frag1_stat[7] + frag2_stat[7] + frag3_stat[7] + frag4_stat[7]
    Bonds               = VLa_stats[8] + VHa_stats[8] + VLb_stats[8] + VHb_stats[8] + frag1_stat[8] + frag2_stat[8] + frag3_stat[8] + frag4_stat[8]
    Bonds               = [x for x in Bonds if x != []]
    disulphide_bridges        = [VHa_stats[11] + VHb_stats[11]] + [VHa_stats[12] + VHb_stats[12]] + [VHa_stats[13] + VHb_stats[13]] + [VHa_stats[14] + VHb_stats[14]]+ [VHa_stats[15] + VHb_stats[15]] + completed_disulphidebridges
    disulphide_bridges        = [x for x in disulphide_bridges if x != []]
    Label_Text          = VLa_stats[9] + VHa_stats[9] + VLb_stats[9] + VHb_stats[9] + frag1_stat[9] + frag2_stat[9] + frag3_stat[9] + frag4_stat[9]
    Label_spot          = VLa_stats[10] +VHa_stats[10] +VLb_stats[10] +VHb_stats[10]+ frag1_stat[10]+ frag2_stat[10]+ frag3_stat[10]+ frag4_stat[10]
    Domain_Text         = VLa_stats[17] + VHa_stats[17] + VLb_stats[17] + VHb_stats[17] + frag1_stat[17] + frag2_stat[17] + frag3_stat[17] + frag4_stat[17]
    Notes               = VHa_stats[18] + VLa_stats[18] + VHb_stats[18] + VLb_stats[18] + frag1_stat[18] + frag2_stat[18] + frag3_stat[18] + frag4_stat[18]
    Notes_positions     = VHa_stats[19] + VLa_stats[19] + VHb_stats[19] + VLb_stats[19] + frag1_stat[19] + frag2_stat[19] + frag3_stat[19] + frag4_stat[19]
    arcs_left           = VHa_stats[20] + VLa_stats[20] + VHb_stats[20] + VLb_stats[20] + frag1_stat[20] + frag2_stat[20] + frag3_stat[20] + frag4_stat[20]
    arcs_right          = VHa_stats[21] + VLa_stats[21] + VHb_stats[21] + VLb_stats[21] + frag1_stat[21] + frag2_stat[21] + frag3_stat[21] + frag4_stat[21]
    arcs_left_slant     = VHa_stats[22] + VLa_stats[22] + VHb_stats[22] + VLb_stats[22] + frag1_stat[22] + frag2_stat[22] + frag3_stat[22] + frag4_stat[22]
    arcs_right_slant    = VHa_stats[23] + VLa_stats[23] + VHb_stats[23] + VLb_stats[23] + frag1_stat[23] + frag2_stat[23] + frag3_stat[23] + frag4_stat[23]
    ADCs                = VHa_stats[24] + VLa_stats[24] + VHb_stats[24] + VLb_stats[24] + frag1_stat[24] + frag2_stat[24] + frag3_stat[24] + frag4_stat[24]
    print(disulphide_bridges)
    print(Light_Domains_a, Light_Domains_b)
    return(Heavy_Domains_a,Light_Domains_a,Heavy_Domains_b,Light_Domains_b,Heavy_Domains_c,Light_Domains_c,Heavy_Domains_d,Light_Domains_d,Bonds,disulphide_bridges,Label_Text,Label_spot,Domain_Text,Notes,Notes_positions,arcs_left,arcs_right,arcs_left_slant,arcs_right_slant,ADCs)

def render(chains_list,canvas):
    canvas.delete("all")

    Heavy_Domains_a    = chains_list[0]
    Light_Domains_a    = chains_list[1]
    Heavy_Domains_b    = chains_list[2]
    Light_Domains_b    = chains_list[3]
    Heavy_Domains_c    = chains_list[4]
    Light_Domains_c    = chains_list[5]
    Heavy_Domains_d    = chains_list[6]
    Light_Domains_d    = chains_list[7]
    Bonds              = chains_list[8]
    disulphide_bridges = chains_list[9]
    Label_Text         = chains_list[10]
    Label_positions    = chains_list[11]
    Domain_Text        = chains_list[12]
    Notes              = chains_list[13]
    Note_positions     = chains_list[14]
    arcs_left          = chains_list[15]
    arcs_right         = chains_list[16]
    arcs_left_slant    = chains_list[17]
    arcs_right_slant   = chains_list[18]
    ADCs               = chains_list[19]



#disulphide_bridge
    if disulphide_bridges != []:
        for i in range(len(disulphide_bridges)):
            canvas.create_line(disulphide_bridges[i], fill='#FF4040', width = 2)

#Bonds
    for i in range(len(Bonds)):
        canvas.create_line(Bonds[i], fill='#000000', width = 2)
    if arcs_left!=[]:
        for i in range(len(arcs_left)):
            canvas.create_arc(arcs_left[i], start=90, extent=180, style=tk.ARC, fill='#000000', width = 2)
    if arcs_right!=[]:
        for i in range(len(arcs_right)):
            canvas.create_arc(arcs_right[i], start=270, extent=180, style=tk.ARC, fill='#000000', width = 2)
    if arcs_left_slant != []:
        for i in range(len(arcs_left_slant)):
            canvas.create_arc(arcs_left_slant[i], start=150, extent=120, style=tk.ARC,width=2)
    if arcs_right_slant != []:
        for i in range(len(arcs_right_slant)):
            canvas.create_arc(arcs_right_slant[i], start=270, extent=120, style=tk.ARC,width=2)
#ADCs
    if ADCs != []:
        for i in range(len(ADCs)):
            canvas.create_polygon(ADCs[i], outline='#000000',fill='#68C1C1', width=2)
#A domains
    for i in range(len(Heavy_Domains_a)):
        canvas.create_polygon(Heavy_Domains_a[i], outline='#000000',fill='#007ECB', width=2)
    for i in range(len(Light_Domains_a)):
        canvas.create_polygon(Light_Domains_a[i], outline='#000000',fill='#73CAFF', width=2)

#B domains
    for i in range(len(Heavy_Domains_b)):
        canvas.create_polygon(Heavy_Domains_b[i], outline='#000000',fill='#FF43EE', width=2)
    for i in range(len(Light_Domains_b)):
        canvas.create_polygon(Light_Domains_b[i], outline='#000000',fill='#F9D3F5', width=2)

#C domains
    for i in range(len(Heavy_Domains_c)):
        canvas.create_polygon(Heavy_Domains_c[i], outline='#000000',fill='#0BD05A', width=2)
    for i in range(len(Light_Domains_c)):
        canvas.create_polygon(Light_Domains_c[i], outline='#000000',fill='#B9FAD3', width=2)
#
#D domains
    for i in range(len(Heavy_Domains_d)):
        canvas.create_polygon(Heavy_Domains_d[i], outline='#000000',fill='#D9DE4A', width=2)
    for i in range(len(Light_Domains_d)):
        canvas.create_polygon(Light_Domains_d[i], outline='#000000',fill='#E2F562', width=2)

    for i in range(len(Label_positions)):
        canvas.create_text(Label_positions[i], text=Domain_Text[i])
    for i in range(len(Note_positions)):
        canvas.create_text(Note_positions[i],text=Notes[i])

    #canvas.pack(fill=BOTH, expand=1)


def render_pipeline(canvas):
    entry=textBox.get("1.0","end-1c")
    split_chains = Get_dictionaries(entry)
    coordinates  = Check_interactions(split_chains)
    render(coordinates, canvas)

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




###############Main programme#######################

root = tk.Tk()

canvas = tk.Canvas(root, height=HEIGHT, width=WIDTH, bg='#E7E0E6')
canvas.pack()

frame = tk.Frame(root, bg = '#80c1ff',bd=5)
frame.place(relwidth=1,relheight=1)

###Input box
textBox = tk.Text(frame, font=40)
textBox.place(relx=0.01,rely = 0.05, relwidth=0.4,relheight=0.3)

###Option box
frame2 = tk.Frame(frame, bg = '#D3D3D3')
frame2.place(relx=0.01, rely = 0.45,relheight = 0.35, relwidth = 0.4)

status_label = tk.Label(root, text='test', bd=1)
status_label.place(rely = 0.98, relheight = 0.02, relwidth = 1)

##Big button
button = tk.Button(frame, text = "Get Structure", bg = "grey", font=40, command=lambda: render_pipeline(lower_canvas))
button.place(relx=0.1,rely=0.85,relheight=0.1, relwidth=0.2)
button.bind("<Enter>", button_hover)
button.bind("<Leave>", button_hover_leave)

###Results canvas
lower_frame = tk.Frame(root, bg = '#80c1ff', bd=10)
lower_frame.place(relx=0.45, rely=0.05, relwidth=0.55,relheight=0.9)

lower_canvas = tk.Canvas(lower_frame)
lower_canvas.place(relheight=1,relwidth=1)




root.mainloop()
