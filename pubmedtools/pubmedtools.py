#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
PubmedTools (pubmedtools) package provides functions for searching and
retrieving articles from the PubMed database using Biopython and NCBI Entrez
Direct. This is not an official NCBI library and has no direct affiliation
with the organization.

Functions:
- pubmed_search_biopython: Searches the PubMed database using a given term and
retrieves article information such as abstract, title, publication date,
authors, MeSH terms, and other terms related to each article. This function
uses the Bio.Entrez package from Biopython. Please note that this function
has a limitation of retrieving a maximum of 10,000 results.
- pubmed_search_edirect: Searches the PubMed database using a given term and
retrieves article information such as abstract, title, publication date,
authors, MeSH terms, and other terms related to each article. This function
uses the official NCBI Entrez Direct tool. This function works with Linux and
Windows systems using WSL (Windows Subsystem for Linux).
- prepare_edirect_folder: Function to prepare the edirect folder for
pubmed_search_edirect. It checks if the edirect folder exists and contains the
necessary files. If not, it downloads and extracts the required files into the
edirect folder.

Author: Diogo de J. S. Machado
Date: 13/07/2023
"""
import re
from Bio import Entrez, Medline
import time
import pandas as pd
from io import StringIO
import platform
import os
import subprocess
import gzip
import zipfile
import urllib.request
import shutil
    
def pubmed_search_biopython(term, email='', api_key=None, batch_size=1000, verbose=True):
    """
    Searches the PubMed database using a given term and retrieves the abstract, title, publication date,
    authors, MeSH terms, and other terms related to each article. This function use the Bio.Entrez package
    from Biopython.
    The search is limited to 10,000 results.

    Parameters
    ----------
    term : str
        The search term to be used in the query.
    email : str, optional
        Email address to be used in case the Entrez server needs to contact you.
    api_key : str, optional
        The NCBI API key. If not provided, the search will be performed without the API key.
    batch_size : int, optional
        Number of articles to be downloaded per iteration. Default is 1000.
    verbose : bool, optional
        Whether to print progress messages. Default is True.

    Returns
    -------
    result : pandas.DataFrame
        A DataFrame with columns 'pmid', 'ti', 'ab', 'fau', 'dp', 'mh', and 'ot'.
        Each row contains information related to a single article retrieved from the search term query.

    Raises
    ------
    Exception
        If the search returns more than 10,000 results, which is the limit of this function.
        In this case, the user should use the `pubmed_search_edirect` function.
    """

    Entrez.email = email
    
    if api_key:
        Entrez.api_key = api_key

    if verbose:
        print(term)

    # ESearch to get webenv and query key
    h = Entrez.esearch(db='pubmed', term=term, usehistory="y")
    search_results = Entrez.read(h)
    count = int(search_results['Count'])
    
    if count > 10000:
        raise Exception(f"This search has {count} results. "+
                        "The pubmed_search_entrez function does not support more than 10000 results. "+
                        "Use pubmed_search_edirect.")
    
    r=[]

    for start in range(0, count, batch_size):
        if verbose:
            print('Downloading %i-%i/%i'%(start, min([start+batch_size,count]), count))

        # EFetch to retrieve data for a given set of PMIDs
        fetch_handle = Entrez.efetch(
            db="pubmed",
            rettype="medline",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )
        
        # Parsing the downloaded data in medline format
        records_medline = Medline.parse(fetch_handle)
        
        for record in records_medline:
            # Extracting the relevant fields for each article
            ab = record.get('AB', '')
            ab = re.sub('\s+',' ',ab)
            ti = record.get('TI', '')
            ti = re.sub('\s+',' ',ti)
            pmid = record.get('PMID', '')
            fau = record.get('FAU', '')
            dp = record.get('DP', '')
            mh = record.get('MH', '')
            ot = record.get('OT', '')
            
            # Creating a dictionary with the relevant fields for each article and adding it to the list
            row = {
            'pmid': int(pmid),
            'ti': ti,
            'ab': ab,                 
            'fau': fau,
            'dp': dp,
            'mh': mh,
            'ot': ot,
            }
           
            r.append(row)
            
        # If not all results have been downloaded, wait for 3 seconds before the next download,
        # to avoid overloading the NCBI server
        if (start+batch_size) < count:
            time.sleep(3)    
    
    print('Done!')
    
    result = pd.DataFrame(data=r)
    return result

def pubmed_search_edirect(query, api_key=None):
    """
    Searches the PubMed database using a given term and retrieves the abstract, title, publication date,
    authors, MeSH terms, and other terms related to each article. This function use the official NCBI
    Entrez Direct tool.

    Parameters
    ----------
    query : str
        The query to be searched in PubMed.
    api_key : str, optional
        The NCBI API key. If not provided, the search will be performed without the API key.

    Returns
    -------
    result : pandas.DataFrame
        A pandas DataFrame containing the search results.

    Notes
    -----
    This function works with Linux and Windows systems using WSL (Windows Subsystem for Linux).

    Raises
    ------
    Exception
        If the operating system is not recognized.
    """

    edirect_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'edirect')

    if not os.path.isdir(edirect_path):
        raise ImportError('The "edirect" folder was not found. Please run the "prepare_edirect_folder()" function to create it.')

    print('Downloading data from PubMed...')

    # Check the operating system where the function is being executed
    if platform.system() == "Linux":
        # Define the command to execute the PubMed search
        cmd = 'esearch -db pubmed -query "%s" | efetch -format medline' % query

        # Define the PATH environment variable to include the "edirect" directory
        env = {'PATH': f'{edirect_path}:{os.environ["PATH"]}'}

        # Check if an NCBI API key has been provided
        if api_key:
            # If an API key has been provided, add it to the command
            env['NCBI_API_KEY'] = api_key

        print(cmd)

        # Execute the command and capture the result as a string
        medline_str = subprocess.check_output(cmd, env=env, shell=True).decode('utf-8')

    elif platform.system() == "Windows":
        # Convert the "edirect" directory path to WSL format
        edirect_path = re.sub('\\\\', '/', re.sub('^(\w)\:\\\\',
                                                 lambda match: '/mnt/' + match.group(1).lower() + '/', edirect_path))

        # Define the command to execute the PubMed search
        cmd = 'wsl {}/esearch -db pubmed -query "{}" | wsl {}/efetch -format medline'.format(edirect_path, query, edirect_path)

        # Print the command to be executed
        print(cmd)

        # Check if an NCBI API key has been provided
        # If yes, add the NCBI_API_KEY environment variable to the command
        if api_key:
            cmd = re.sub('wsl', f'wsl export NCBI_API_KEY={api_key};', cmd)

        # Execute the command and capture the result as a string
        medline_str = subprocess.check_output(cmd, shell=True).decode('utf-8')

    else:
        # Raise an error if the operating system is not recognized
        raise Exception('Unsupported operating system: {}'.format(platform.system()))

    # Extract relevant data from the search result
    medline_record = Medline.parse(StringIO(medline_str))
    result = []
    print('Extracting data...')
    for record in medline_record:
        ab = record.get('AB', '')
        ti = record.get('TI', '')
        pmid = record.get('PMID', '')
        aut = record.get('FAU', '')
        dp = record.get('DP', '')
        mh = record.get('MH', '')
        ot = record.get('OT', '')

        # Store the relevant data in a dictionary
        row = {
            'pmid': int(pmid),
            'ti': ti,
            'ab': ab,
            'fau': aut,
            'dp': dp,
            'mh': mh,
            'ot': ot,
        }
        result.append(row)

    # Create a DataFrame from the dictionary
    print('Writing file...')
    result = pd.DataFrame(data=result)

    print('Done!')

    # Return the final result
    return result

def prepare_edirect_folder():
    """
    Function to prepare the edirect folder for pubmed_search_edirect.

    Checks if the edirect folder exists and contains the necessary files.
    If not, it downloads and extracts the required files into the edirect folder.
    """
    # Get the module directory and create the edirect path
    module_dir = os.path.dirname(os.path.abspath(__file__))
    edirect_path = os.path.join(module_dir, "edirect")

    # Check if the edirect folder exists, if not, create it
    if not os.path.exists(edirect_path):
        os.makedirs(edirect_path)

    # Check if the esearch file is present in the edirect folder
    if not os.path.exists(os.path.join(edirect_path, "esearch")):
        print("Downloading and extracting edirect...")

        # Download the edirect ZIP file
        edirect_url = "https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip"
        edirect_zip_path = os.path.join(edirect_path, "edirect.zip")
        urllib.request.urlretrieve(edirect_url, edirect_zip_path)

        # Extract the contents of the edirect ZIP file
        with zipfile.ZipFile(edirect_zip_path, "r") as zip_ref:
            zip_ref.extractall(edirect_path)

        # Get the extracted edirect directory
        edirect_extracted_dir = os.path.join(edirect_path, "edirect")

        # Get the list of files in the extracted directory
        files = os.listdir(edirect_extracted_dir)

        # Move each file to the edirect folder
        for f in files:
            origin_path = os.path.join(edirect_extracted_dir, f)
            dest_path = os.path.join(edirect_path, f)
            shutil.move(origin_path, dest_path)

        # Remove the edirect ZIP file and the extracted directory
        os.remove(edirect_zip_path)
        os.rmdir(edirect_extracted_dir)

        # Download xtract.Linux.gz
        xtract_url = "https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz"
        xtract_gz_path = os.path.join(edirect_path, "xtract.Linux.gz")
        urllib.request.urlretrieve(xtract_url, xtract_gz_path)

        # Extract xtract.Linux.gz
        with gzip.open(xtract_gz_path, "rb") as f_in:
            with open(os.path.join(edirect_path, "xtract"), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Remove xtract.Linux.gz
        os.remove(xtract_gz_path)

        print("EDirect ready!")
    else:
        print("EDirect already ready!")
