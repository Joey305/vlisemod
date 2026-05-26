VIRAL_DATABASE

Project Overview
The VIRAL_DATABASE project is an initiative by the Schurer Lab to create a comprehensive repository of viral protein structures from various databases like RCSB PDB. The goal is to facilitate the research and analysis of viral proteins, specifically focusing on their interactions with various ligands, to better understand their mechanisms and potentially guide the development of antiviral drugs.

Deployment defaults
The public/default V-LiSEMOD configuration keeps `PROTACability` enabled in the main navigation, keeps `PROTAC Builder` available as a standalone external link, hides the `Drug GPT` navigation item, and disables Drug GPT plus local LLM loading by default. No model weights should load or download unless you explicitly opt in with environment variables.

Environment toggles
Set `SHOW_DRUG_GPT_NAV=1` to show the Drug GPT nav button.
Set `ENABLE_DRUG_GPT=1` to enable the Drug GPT routes and blueprint.
Set `ENABLE_LOCAL_LLM=1` to allow local LLM configuration and model loading.

Objectives
Data Collection: Automate the retrieval of viral protein structures from RCSB PDB using Python scripts.
Data Curation: Identify and segregate proteins based on the presence or absence of specific ligands.
Data Analysis: Perform structural and surface area analyses to identify potential targets for drug binding.
Database Creation: Compile the curated data into a structured database that can be easily accessed and utilized by researchers.
Repository Contents
1A_PDBpull_1A.py: Script to pull PDB files based on specific viral names.
1B_PDBpull_1B.py: Extended script to pull PDB files for multiple variations of viral names.
2_PDBSORT_2.py: Script to sort PDB files into those with and without specified ligands.
3_SASA_3.py: Script to analyze the solvent-accessible surface area of the proteins.
4_Database_4.py: Script to integrate all data into a MySQL database for easy querying.
reassemble.py: Script to reconstruct the viral database from preprocessed CSV files.
Database_DATA: Directory where all downloaded and processed PDB files are stored.
Installation
Prerequisites
Python 3.6 or higher
Required dependencies:
bash
Copy code
pip install -r requirements.txt
This will install all the necessary dependencies to run the scripts.

Reconstructing the Database
After cloning the repository and installing the required dependencies, follow these steps to reconstruct the viral database:

Run the Reassemble Script:

This script will reconstruct the database from the preprocessed CSV files.
Use the following command to run it:
bash
Copy code
python reassemble.py
The script will populate the database with all the necessary tables and data.
Launch the Web Application:

After the database is successfully reconstructed, you can launch the Flask web application by running:
bash
Copy code
flask run
This will start the local server and make the web app accessible in your browser.
Usage
To use the individual scripts, navigate to the project directory and run:

bash
Copy code
python <script_name>.py

Replace <script_name> with the name of the script you want to run. Make sure to configure any necessary parameters such as API keys or database credentials in the scripts or in a configuration file.

Contributing
We welcome contributions from the scientific and developer community. If you have suggestions or improvements, please fork the repository and submit a pull request.

License
This project is open-sourced under the MIT license.

Contact
For questions or collaboration requests, please contact the Schurer Lab at Sschurer@med.miami.edu.

Acknowledgments
We thank the contributors and supporters of the VIRAL_DATABASE project for their insights and dedication to advancing viral protein research.
