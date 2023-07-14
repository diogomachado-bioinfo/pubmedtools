# PubmedTools (pubmedtools)
PubmedTools (pubmedtools) package provides functions for searching and
retrieving articles from the PubMed database using Biopython and NCBI Entrez
Direct. This is not an official NCBI library and has no direct affiliation
with the organization.

---

---
# Features

- `pubmedtoos.search.biopython_search`: Searches the PubMed database using a
                                               Biopython Entrez module (Bio.Entrez).
- `pubmedtoos.search.edirect_search`: Searches the PubMed database using the
                                             official Entrez Direct tool.
- `pubmedtoos.prepenv.edirect_folder`: Prepares the Entrez Direct folder for use
                                       with the edirect_search function.

---

---
# Installation
You can install PubmedTools using pip:

    pip install pubmedtools

---

---
# Functions

## search

---
### `pubmedtools.search.biopython_search`
Searches the PubMed database using a given term and retrieves the abstract,
title, publication date, authors, MeSH terms, and other terms related to each
article. This function use the Bio.Entrez module from Biopython. The search is
limited to 10,000 results.

**Parameters**

- `term` : str
    - The search term to be used in the query.
- `email` : str, optional
    - Email address to be used in case the Entrez server needs to contact you.
- `api_key` : str, optional
    - API key to access the Entrez server.
- `batch_size` : int, optional
    - Number of articles to be downloaded per iteration. Default is 1000.
- `verbose` : bool, optional
    - Whether to print progress messages. Default is True.

**Returns**

- pandas.DataFrame
    - A DataFrame with columns 'pmid', 'ti', 'ab', 'fau', 'dp', 'mh', and 'ot'.
    - Each row contains information related to a single article retrieved from
    the search term query.

**Raises**

Exception
    - If the search returns more than 10,000 results, which is the limit of
    this function.
    In this case, the user should use the `pubmedtools.search.edirect_search`
    function.

### `pubmedtools.search.edirect_search`
Searches the PubMed database using a given term and retrieves the abstract,
title, publication date, authors, MeSH terms, and other terms related to each
article. This function use the official NCBI Entrez Direct tool.

**Parameters**

- `query` : str
    - The query to be searched in PubMed.
- `api_key` : str, optional
    - The NCBI API key. If not provided, the search will be performed without
    the API key.

**Returns**

- pandas.DataFrame
    - A pandas DataFrame containing the search results.

**Notes**

- This function works with Linux and Windows systems using WSL (Windows
Subsystem for Linux).

**Raises**

- Exception
    - If the operating system is not recognized.

---
## `prepenv`

---
### `pubmedtools.prepenv.edirect_folder`
Function to prepare the edirect folder for pubmed_search_edirect. Checks in
pubmedtools package path if the edirect folder exists and contains the
necessary files. If not, it downloads and extracts the required files.