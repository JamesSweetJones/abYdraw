abYdraw

This is a programme designed to use our group's Antibody Markup Language (AbML) for describing bispecific antibody (BsAb) formats by either inputting an AbML descriptor string of a BsAb or by drawing a BsAb and outputting the its descriptor string. It is written in Python 3 and using standard packages TKinter to build the graphical interface in order to make it as accessible as possible.

Contents
1. Installing and Executing
2. Interface
3. AbML
4. Inputting AbML
5. Obtaining AbML
6. Formats Library
7. Saving and exporting
8. Settings



1. Installing and Executing


abYdraw may be downloaded and ran as an executable Python script.


2. Interface


The programme interface includes six points of reference, four of which in a column on the left hand side and two more on the right hand side.
Starting with the left hand column, the first is the Domain palette which has buttons necessary for drawing antibody domains, secondly a library of commonly used bispecific antibody AbML expressions, thirdly the input box for AbML expressions and a buttonpad that will render antibody schematics or output AbML to the textbox.
On the right hand side, the most prominent feature is the canvas for drawing and rendering antibody schematics and underneath there are two buttons which are involved in exporting the schematic.

3. AbML


Our language was derived from existing macromolecule descriptor languages but we have compensated for their limitations and made AbML simple whilst conveying as much useful information as possible. Strings are split into chains, which are then split into domains. Each domain type has its own symbol and each domain unit also carries additional information including: modification types; the specificity of the variable region (if applicable); a number label assigned to the domain and the number label assigned to the domain it interacts with; the number of disulphide bonds between the two interacting domains and comments outlining additional information not covered by the language of types: TYPE; NOTE; MOD; ANTI and LENGTH. Full descriptions of AbML can be found on the language guide sheet included in the Repository.

4. Inputting AbML


AbML descriptor strings my be inputted in the entry box or opened in the "File>Open" menu and then by clicking "Get Structure", a schematic of that antibody will render in the canvas. Schematics are drawn in colour-coded fashion depending on any specificities given in the descriptor chain. Domains are connected by different kinds of linkers which are also colour-coded depending on their type. Any comments given in the descriptor string are also displayed beside the schematic. Labels on the schematic may be toggled on and off using the Labels key on the Domain Palette.

5. Obtaining AbML


To draw an antibody, you must insert domains onto the canvas and arrange them so the programme recognises it as an antibody. Your tools for adding domains to the canvas are in the Domain Palette which contains all of the domain types, modifications, specificities and comment types as described in the AbML guidesheet as well as some options. Selecting a button on the palette will cause it to flash red to indicate it has been switched on. Only one domain or connector type may be switched on at a time, but you may choose any combination of modifications and specificities to accompany your selection. Specificity types are only applicable to VH or VL domains but when drawing other domains, a default "a" value is set. Once selected you will notice the cursor will change from arrow to "+" sign. This means you can left-click to insert your chosen domain type onto the canvas at the location you have clicked. If no domain types are selected but a modification or specificity are, then by clicking on a domain on the canvas, you may replace its current specificity and modification to those you have selected.
Domains may be connected with the connector options in the first column of the palette. When a connector type is selected it will become highlighted in red and the user must click and drag the bond from its starting position to its end, making sure each end is inside the boundaries two domains it links. Bonds are unidirectional and start from N-terminus to C-terminus.
Comments may be added by selecting a comment type, which will highlight the button just pressed and the "Comment" button and then inputting the comment into the entry box beneath the palette. Comments may then be drawn on the domains they are applicable to. To disable further commenting, ensure the comment type and "Comment" buttons are no longer highlighted.
If no domain types or modifications are selected then domains, linkers and comments may be rearranged by clicking and dragging them. You must ensure that interacting domains are positioned close together at roughly the same level and that VH and VL domains face each other to complete their antigen binding domains. Their orientations can be changed by right-clicking the relevant domains. Furthermore features may be deleted by selecting the "Delete" button on the palette and then selecting what to delete. To remove all features, click "Clear All" and the canvas will be made blank.
Once domains and connectors are arranged, click the "Get AbML" button to generate the AbML descriptor string for this sequence. Once this has been generated it will appear in the input box. You may then click "Get structure" again to re-render the schematic with abYdraw. Alternatively, the "Tidy" button performs both steps of this operation. Once rendered, an image may be altered by adding or removing domains. By clicking "Get Sequence" or "Tidy", you will obtain a new expression for rendering.

6. AbML Formats Library


To assist users, the programme has a library of BsAb formats available which can be scrolled through and selected. This will give the schematic and AbML expression for this format that can be used as a starting point to make new expressions and schematics that are relevant to the user.

7. Exporting and Saving


Exporting the canvas image as .eps file can be done by "Export EPS" and using the file directory to save the image. Alternatively the AbML may be saved as a text file by clicking the "File>Save" option in the menu.
Finally you may export a Template File from the expression in the entry box which is a format of notating important BsAb residues. The programme cannot locate these residues but it can identify the features that are in the BsAb that users may want to include in the Template Files.

8. Settings

Users may change aspects about the rendering and pairing of their schematic. By opening the "File>Settings" window, a menu with two tabs will appear. The first tab has settings regarding pairing sensitivity of drawn schematics, bond thickness, directional arrows, Hinge and Linker labels which can be set by using the appropriate sliders. For pairing sensitivity the scale is 0-80 pixels and bond width are between 1-5 pixels. Other binary settings are on sliders 0-1. To update the settings, users must press the update button and re-render their schematic to see their new schematic. The second tab is the  colour-coding  menu which with a list of domain types. When a domain type is selected the current colour of assigned to that domain will appear. "Change colour" allows users to assign a new colours to that domain type using the colour palette of the operating system, but these may be reverted by "Revert colour" or "Revert all colours".


Software authors: James Sweet-Jones and Andrew Martin (Darwin Building, University College London, Gower Street, London)
