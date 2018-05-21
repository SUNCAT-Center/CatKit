## Using the cathub cli

Run `cathub`, like so

    cathub --help

or with any of its sub-commands, like so

    cathub make_folders_template --help

## Examples

To create an .json input file

    cathub make_folders_template project1.json --create-template

To create a folder structures from a .json input file

    cathub make_folders_template project1.json

Querying the Catalysis Hub database:

    cathub reactions -q reactants=CO -q chemicalComposition=~Pt

    cathub publications -q title=~Evolution -q year=2017

Reading folders into sqlite3 db file:

    cathub folder2db <foldername>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
