# abYdraw

This is a programme designed to use our group's Antibody Markup Language (AbML) for describing bispecific antibody (BsAb) formats by either inputting an AbML descriptor string of a BsAb or by drawing a BsAb and outputting the its descriptor string. It is written in Python 3 and with both command-line and graphical interfaces built by using standard packages TKinter to in order to make it as accessible as possible.

### Contents

1. AbML
2. Installing and Executing
3. Interface
4. Inputting AbML
5. Obtaining AbML
6. Formats Library
7. Saving and exporting
8. Settings

Reference:

Software authors: James Sweet-Jones and Andrew Martin (Darwin Building, University College London, Gower Street, London)


### 1. AbML

![AbML Guide Sheet](https://github.com/JamesSweetJones/abYdraw/blob/main/AbML_Guide_Sheet.png)
Our language was derived from existing macromolecule descriptor languages but we have compensated for their limitations and made AbML simple whilst conveying as much useful information as possible. Strings are split into chains, which are then split into domains. Each domain type has its own symbol and each domain unit also carries additional information including: modification types; the specificity of the variable region (if applicable); a number label assigned to the domain and the number label assigned to the domain it interacts with; the number of disulphide bonds between the two interacting domains and comments outlining additional information not covered by the language of types: `TYPE; NOTE; MOD; ANTI and LENGTH`. TYPE and MOD are limited to reserved keywords in block capitals while other comments are not restricted and written in lower case. Full descriptions of AbML can be found on the language guide sheet included in the Repository.

AbML expression for a standard IgG molecule would be:
`VH.a(1:6)-CH1(2:7){1}-H(3:10){2}-CH2(4:11) - CH3(5:12) | VL.a(6:1)-CL(7:2){1} | VH.a(8:13)- CH1(9:14){1}-H(10:3){2}-CH2(11:4)-CH3(12:5) | VL.a(13:8)-CL(14:9){1}`

Full documentation of AbML is available in the [supplementary material](https://github.com/JamesSweetJones/abYdraw/blob/main/paper/SupplementaryMaterial1.md) of the publication 

### 2. Installing and Executing

##### Command Line Interface

abYdraw has a Command Line based interfaced which may be accessed by executing the script with suitable arguements.

usage:
`python abYdraw [options]`

optional arguments:
 
-  -f / --file      the path to plaintext file with AbML expression
-  -i / --input     string of AbML input
-  -o / --output    string of image output name      (default "abYdraw_export")
-  -s / --show      Show image window            0-1 (default 0)
-  -e / --format    Specify image format         "eps", "png", "jepg" (default eps)
-  -i / --image     Save image file              0-1 (default 1)
-  -t / --template  Save template file           0-1 (default 0)
-  -l / --labels    Toggle domain labels         0-1 (default 1)
-  -j / --hinge     Toggle hinge labels          0-1 (default 0)
-  -k / --linker    Toggle linker labels         0-1 (default 0)
-  -a / --arrows    Toggle bond direction arrows 0-1 (default 0)
-  -b / --thickness Set bond thickness           1-5 (default 2)
-  -h / --help      Show this help message and exit

either -f or -i arguements are required to run the script and other options may be toggled between 0 (off) and 1 (on).

##### Graphical Interface

abYdraw may be downloaded and the GUI may be opened by executing the Python script without arguements.

`python <path/to/file>/abYdraw.py`

##### Compiling abYdraw

However, it may also be compiled into an executable file using py2app and the setup.py file provided. which may then be moved into the user's operating system's applications folder for convenient startup.
install a virtual environment 

`pip install virtualenv`

`virtualenv venv --system-site-packages`

Then open a virtual environment and install py2app

`source venv/bin/activate`

`pip install -U py2app`

Finally, run py2app with the setup.py file

`python <path/to/file>/setup.py py2app -A`

This command will generate ./build ./dist folders. Your compiled app will be in the ./dist folder.

Precompiled versions of abYdraw are available for Windows and MacOS on our research group's [website](http://www.bioinf.org.uk/software/abydraw/)

### 3. Graphical Interface

![abYdraw Interface](https://github.com/JamesSweetJones/abYdraw/blob/main/Interface_new.png)

The programme interface includes six points of reference, four of which in a column on the left hand side and two more on the right hand side. Starting with the left hand column, the first is the Domain palette **(A)** which has buttons necessary for drawing antibody domains, secondly a library of commonly used bispecific antibody AbML expressions **(B)**, thirdly the input box for AbML expressions **(C)** and a buttonpad that will render antibody schematics or output AbML to the textbox **(D)**. On the right hand side, the most prominent feature is the canvas for drawing and rendering antibody schematics **(E)** and underneath there are two buttons which are involved in exporting the schematic **(F)**.

### 4. Inputting AbML

AbML descriptor strings my be inputted in the entry box or opened in the `File>Open` menu and then by clicking `Get Structure`, a schematic of that antibody will render in the canvas. Schematics are drawn in colour-coded fashion depending on any specificities given in the descriptor chain. Domains are connected by different kinds of linkers which are also colour-coded depending on their type. Any comments given in the descriptor string are also displayed beside the schematic. Labels on the schematic may be toggled on and off using the `Labels` key on the Domain Palette.

### 5. Obtaining AbML

##### Drawing BsAb Domains
To draw an antibody, you must insert domains onto the canvas and arrange them so the programme recognises it as an antibody. Tools for adding domains to the canvas are in the Domain Palette which contains all of the domain types, modifications, specificities and comment types as described in the AbML guidesheet as well as some options. Selecting a button on the palette will cause it to flash red to indicate it has been switched on. Only one domain or connector type may be switched on at a time, but you may choose any combination of modifications and specificities to accompany your selection. Specificity types are only applicable to variable domains where if no specificity is selected, it will be rendered as a default `a` value. Mulitple specificities may be selected by right-clicking the first specificity and left-clicking subsequent specificities. Once a domain is selected you will notice the cursor will change from arrow to "+" sign. This means you can left-click to insert your chosen domain type onto the canvas at the location you have clicked. 

##### Connecting Domains
Domains may be connected with the connector options in the first column of the palette. When a connector type is selected it will become highlighted in red and the user must click and drag the bond from its starting position to its end, making sure each end is inside the boundaries two domains it links. Bonds are unidirectional and start from N-terminus to C-terminus.

##### Domain Comments
Comments may be added by selecting a comment type, which will highlight the button just pressed and the `Comment` button and then inputting the comment into the entry box beneath the palette, However, `TYPE` or `MOD` comments must be selected from the appropriate drop-down boolean lists because these are reserved values in AbML. Comments may then be drawn on the domains they are applicable to. To disable commenting, ensure the comment type and `Comment` buttons are no longer highlighted. 

##### Editing Drawings
If no domain types or modifications are selected then but a modification or specificity are, then by clicking on a domain on the canvas, you may replace its current specificity and modification to those that are currently selected. To remove a modification, specificity or comment from a domain then select the feature you wish to remove from Pallete and click a domain which already has that feature. The feature will then be removed from the Domain and using `Tidy` will remove the features from the AbML string and the schematic.

Additionally if no Domain Pallete buttons are selected, domains, linkers and comments may be rearranged by clicking and dragging them when the curseor is in the "Fleur" configuration (a cross with arrows at all points). You must ensure that interacting domains are positioned close together at roughly the same level and that VH and VL domains face each other to complete their antigen binding domains. If the new position of the domains no longer links any connectors origianlly inside those domains, these original connectors are deleted. Their orientations can be changed by right-clicking the relevant domains. 

Furthermore features may be deleted by selecting the "Delete" button on the palette and then selecting what to delete. To remove all features, click `Clear All` and the canvas will be made blank. 

##### Get AbML

Once domains and connectors are arranged, click the `Get AbML` button to generate the AbML descriptor string for this sequence. Once this has been generated it will appear in the input box. You may then click `Get Structure` again to re-render the schematic with abYdraw. Alternatively, the `Tidy` button performs both steps of this operation. Once rendered, an image may be altered by editing, adding or removing domains. By clicking `Get AbML` or `Tidy`, you will obtain a new expression for rendering.


### 6. AbML Formats Library

To assist users, the programme has a library of BsAb formats available which can be scrolled through and selected. This will give the schematic and AbML expression for this format that can be used as a starting point to make new expressions and schematics that are relevant to the user.

### 7. Exporting and Saving

Exporting the canvas image as .eps file can be done by `Export Image` and using the file directory to save the image in JPEG, PNG or EPS format. Alternatively the AbML may be saved as a text file by clicking the `File>Save` option in the menu. Users may also export a Template File from the expression in the entry box which is a format of notating important BsAb residues. The programme cannot locate these residues but it can identify the features that are in the BsAb that users may want to include in the Template Files.

### 8. Settings

Users may change aspects about the rendering and pairing of their schematic. By opening the `File>Settings` window, a menu with two tabs will appear. The first tab has settings regarding pairing sensitivity of drawn schematics, bond thickness, directional arrows, Hinge and Linker labels which can be set by using the appropriate sliders. For pairing sensitivity the scale is 0-80 pixels and bond width are between 1-5 pixels. Other binary settings are on sliders 0-1. To update the settings, users must press the update button and re-render their schematic to see their new schematic. The second tab is the colour-coding menu which with a list of domain types. When a domain type is selected the current colour of assigned to that domain will appear. `Change colour` allows users to assign a new colours to that domain type using the colour palette of the operating system, but these may be reverted by `Revert colour` or `Revert all colours`.



