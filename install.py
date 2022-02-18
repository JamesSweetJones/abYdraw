#!/usr/bin/python3
#*************************************************************************
#
#   Program:    abYdraw
#   File:       install.py
#   
#   Version:    V1.0
#   Date:       18.02.22
#   Function:   Install the abYdraw software for Linux
#   
#   Copyright:  (c) Prof. Andrew C. R. Martin, UCL, 2022
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   18.02.22  Original   By: ACRM
#
#*************************************************************************
import sys
import subprocess
import os.path

#*************************************************************************
# Main program setup

# Define the list of files we wish to install
files = ['abYdraw.py']
# Define and links we wish to create to those files
links = {'abYdraw.py': 'abydraw'}

#*************************************************************************
def usage_die():
    """Prints a usage message
       18.02.22 Original   By: ACRM
    """
    
    print("\nabYdraw install script")
    print("V1.0 (c) 2022 UCL, Andrew C.R. Martin, James Sweet-Jones")
    print("\nUsage: ./install.py [dest]")
    print("       dest - specify the destination [default ~/bin/]")
    print("\nInstalls the abYdraw script and creates a link such that")
    print("it can be run simply by typing 'abydraw'\n")
    exit(0)

    
#*************************************************************************
def getPackageManager():
    """Identifies the package manager for Linux
       18.02.22 Original   By: ACRM
    """
    if os.path.exists('/usr/bin/dnf') or os.path.exists('/bin/dnf'):
        return('dnf')
    elif os.path.exists('/usr/bin/yum') or os.path.exists('/bin/yum'):
        return('yum')
    elif os.path.exists('/usr/bin/apt-get') or os.path.exists('/bin/apt-get'):
        return('apt-get')
    return('none')

   
#*************************************************************************
def install(rpm_module, apt_module):
    """Install an RPM or APT package - determines which package
       manager is in use and installs the correct version
       18.02.22 Original   By: ACRM
    """
    module = rpm_module

    package_manager = getPackageManager()
    if package_manager == 'apt-get':
        module = apt_module
    elif package_manager == 'none':
        return False
        
    print("Installing module '" + module + "'...")
    exe = "sudo " + package_manager + " -y install " + module
    try:
        retval = subprocess.check_output(exe, shell=True)
        retval = str(retval, 'utf-8')
        print (retval)
    except:
        print ("\nUnable to run package manager")
        print ("Installation of Python module failed")
        exit (1)
    return True


#*************************************************************************
def testandinstall_tkinter():
    """Tests whether the tkinter Python module is installed.
       If not, installs it as a Linux package if there is a
       package manager available. If not, it reports a message
       saying the user will have to do that manually
       18.02.22 Original   By: ACRM
    """
    try:
        import tkinter
    except ModuleNotFoundError:
        if(not install('python3-tkinter', 'python3-tk')):
            print("You will need to install the " + module +
                  "module for Python3 manually")
            exit(1)

            
#*************************************************************************
def copy_files(files, dest):
    """Copies a set of files to the given destination
       18.02.22 Original   By: ACRM
    """
    for file in files:
        exe = "cp " + file + " " + dest
        print (exe)
        try:
            subprocess.check_output(exe, shell=True)
        except:
            if not os.path.exists(dest + "/" + file):
                print ("\nUnable to copy files to directory: " + dest)
                print ("Installation failed")
                exit (1)
    return
    

#*************************************************************************
def make_links(links, dest):
    """Creates a set of symbolic links within the dest dorectory
       18.02.22 Original   By: ACRM
    """
    for src in links:
        lnk = links[src]
        exe = "(cd " + dest + "; ln -sf " + src + " " + lnk + ")"
        print (exe)
        try:
           subprocess.check_output(exe, shell=True)
        except:
            print ("\nUnable to create links in: " + dest)
            print ("Installation failed")
            exit(1)
    return
    

#*************************************************************************
def make_dir(dir):
    """Creates a directory if it doesn't exist
       18.02.22 Original   By: ACRM
    """
    if not os.path.isdir(dir):
        print ("Creating directory " + dir)
        exe = "mkdir -p " + dir
        try:
            subprocess.check_output(exe, shell=True)
        except:
            if not os.path.isdir(dir):
                print ("Unable to create directory: " + dir)
                print ("Installation failed")
                exit (1)
    return

#*************************************************************************
###                        MAIN PROGRAM                                ###
#*************************************************************************
### 18.02.22 Original   By: ACRM

# Find out where we are installing
dest = os.getenv('HOME') + "/bin"
if len(sys.argv) > 1:
    dest = sys.argv[1]

if dest == '-h' or dest == '-help' or dest == '--help':
    usage_die()

testandinstall_tkinter()
make_dir(dest)
copy_files(files, dest)
make_links(links, dest)
