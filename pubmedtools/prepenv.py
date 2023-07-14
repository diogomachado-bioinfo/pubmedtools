#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
This module provides the function to prepare environment to run
pubmedtools.pubmed_search_edirect.
The pubmedtools is not an official NCBI library and has no direct affiliation
with the organization.

Function:
- edirect_folder: Function to prepare the edirect folder for
pubmed_search_edirect. It checks if the edirect folder exists and contains the
necessary files. If not, it downloads and extracts the required files into a
folder in pubmedtools package.

Author: Diogo de J. S. Machado
Date: 14/07/2023
"""
import os
import gzip
import zipfile
import urllib.request
import shutil

def edirect_folder():
    """
    Function to prepare the edirect folder for pubmed_search_edirect. Checks in
    pubmedtools package path if the edirect folder exists and contains the
    necessary files. If not, it downloads and extracts the required files.
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
