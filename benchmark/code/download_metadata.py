import requests
import os
import xml.etree.ElementTree as ET

base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

def downloadGSM4(GSM, GSE):
    gsm_id = GSM
    gse_id = GSE

    # Function to download and save XML file
    def download_xml(id):
        params = {
            'acc': id,
            'targ': 'self',
            'view': 'full',
            'form': 'xml'
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            xml_folder = 'XMLFiles'  # Folder to save XML files
            if not os.path.exists(xml_folder):
                os.makedirs(xml_folder)
            filename = os.path.join(xml_folder, f"{id}.xml")  # Save file in XMLFiles folder
            with open(filename, 'wb') as file:
                file.write(response.content)
            print(f"Saved {id} to {filename}")
            return filename
        else:
            print(f"Failed to retrieve {id}")
            return None

    # Function to safely get text from XML
    def find_element_text(root, path, namespaces, default="Not available"):
        element = root.find(path, namespaces)
        return element.text.strip() if element is not None else default

    # Namespaces for XML parsing
    namespaces = {
        'miniml': 'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'
    }

    # Initialize folder and file paths
    folder_name = 'downloaded'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    output_file_path = os.path.join(folder_name, f'{gsm_id}.txt')

    # Download and parse GSE XML first
    gse_xml = download_xml(gse_id)
    gse_content = ""
    if gse_xml:
        tree = ET.parse(gse_xml)
        root = tree.getroot()

        # Extract GSE information
        gse_title = find_element_text(root, 'miniml:Series/miniml:Title', namespaces)
        gse_summary = find_element_text(root, 'miniml:Series/miniml:Summary', namespaces)
        gse_overall_design = find_element_text(root, 'miniml:Series/miniml:Overall-Design', namespaces)
        sample_type = find_element_text(root, 'miniml:Sample/miniml:Type', namespaces) 

        # Construct GSE information string
        gse_content += f"GSE ID: {gse_id}\nGSE Title: {gse_title}\nExperiment Type: {sample_type}\nGSE Summary: {gse_summary}\nGSE Overall design: {gse_overall_design}\n"


    gsm_xml = download_xml(gsm_id)
    if gsm_xml:
        tree = ET.parse(gsm_xml)
        root = tree.getroot()

        sample_iid = root.find('miniml:Sample', namespaces).get('iid')
        title = find_element_text(root, 'miniml:Sample/miniml:Title', namespaces)
        # ... [other code for extracting GSM information] ...

        sample_type = find_element_text(root, 'miniml:Sample/miniml:Type', namespaces)
        source_name = find_element_text(root, 'miniml:Sample/miniml:Channel/miniml:Source', namespaces)
        organism = find_element_text(root, 'miniml:Sample/miniml:Channel/miniml:Organism', namespaces)

        characteristics_list = root.findall('miniml:Sample/miniml:Channel/miniml:Characteristics', namespaces)
        characteristics = {element.get('tag'): element.text.strip() for element in characteristics_list if element.text}

        treatment_protocol = find_element_text(root, 'miniml:Sample/miniml:Channel/miniml:Treatment-Protocol', namespaces)
        growth_protocol = find_element_text(root, 'miniml:Sample/miniml:Channel/miniml:Growth-Protocol', namespaces)

        # Extract Extracted Molecule and Extraction Protocol
        extracted_molecule = find_element_text(root, 'miniml:Sample/miniml:Channel/miniml:Molecule', namespaces)
        extraction_protocol = find_element_text(root, 'miniml:Sample/miniml:Channel/miniml:Extract-Protocol', namespaces)

        library_strategy = find_element_text(root, 'miniml:Sample/miniml:Library-Strategy', namespaces)
        library_source = find_element_text(root, 'miniml:Sample/miniml:Library-Source', namespaces)
        library_selection = find_element_text(root, 'miniml:Sample/miniml:Library-Selection', namespaces)
        data_processing = find_element_text(root, 'miniml:Sample/miniml:Data-Processing', namespaces)

        # Now construct the string to write to a text file
        output_content = f"GSM ID: {sample_iid}\nTitle: {title}\nSample Type: {sample_type}\nSource Name: {source_name}\nOrganism: {organism}\n"

        # Add Characteristics
        for char_key, char_value in characteristics.items():
            output_content += f"{char_key.capitalize()}: {char_value}\n"

        #add extracted molecule info

        # Add Treatment Protocol and Growth Protocol
        output_content += f"""Treatment Protocol: {treatment_protocol}\nGrowth Protocol: {growth_protocol}
                                            \nExtracted molecule: {extracted_molecule}\nExtraction protocol: {extraction_protocol}
                                            \nLibrary Strategy: {library_strategy}\nLibrary Source: {library_source}\nLibrary Selection: {library_selection}
                                            \nData Processing: {data_processing}\n"""




    with open(output_file_path, 'w', encoding='utf-8') as file:
        combined_content = gse_content + "\n" + output_content
        if "single-cell RNA-Seq" not in combined_content:
            file.write(combined_content)
        else:
            print("The text contains the phrase 'single-cell RNA-Seq' and will not be written to the file.")

    print(f"GSE and GSM information processed and checked for the specific phrase.")



import csv

# Initialize an empty list
mylist4 = []

# getting all the gsm and gse id's to download
with open('goal_to_1M.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile)

    for i, row in enumerate(csvreader):
        if i > 0:
            mylist4.append((row[0], row[1]))



for gsm, gse in mylist4:
    downloadGSM4(gsm, gse)
