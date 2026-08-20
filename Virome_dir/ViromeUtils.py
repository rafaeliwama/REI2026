import pandas as pd


def getNCBIdics(NodeFile, NamesFile):
    """    Get dictionaries from NCBI taxonomy files.
    Parameters:
    NodeFile (str): Path to the nodes.dmp file.
    NamesFile (str): Path to the names.dmp file.
    Returns:
    tuple: A tuple containing four dictionaries:
        - names_dict: Mapping of taxon names to their IDs. {'ID': 'name'}
        - IDnames_dict: Mapping of taxon IDs to their names. {'name': 'ID'}
        - nodesParent_dict: Mapping of taxon IDs to their parent IDs. {'ID': 'parent_ID'}
        - nodesRank_dict: Mapping of taxon IDs to their Linnean rank. {'ID': 'LinneanRank'}
    """
    # Read the names and nodes files into pandas DataFrames

    names_file = pd.read_csv(NamesFile, sep='\t|\t', header=None, engine='python')
    nodes_file = pd.read_csv(NodeFile, sep='\t|\t', header=None, engine='python')
    names_dict = dict(zip(names_file[2], names_file[0]))
    IDnames_dict = dict(zip(names_file[0], names_file[2]))
    nodesParent_dict = dict(zip(nodes_file[0], nodes_file[2]))
    nodesRank_dict = dict(zip(nodes_file[0], nodes_file[4]))

    return names_dict, IDnames_dict, nodesParent_dict, nodesRank_dict

## defining functions

def get_NodepPath(node, NodesParent):
    """This function retrieves the path of parent nodes from a given node to the root node (1).
    It returns a list of nodes in the path."""

    node_path = [node]
    ac_node = node
    while ac_node != 1:
        node_path.append(NodesParent[ac_node])
        ac_node = NodesParent[ac_node]
    return node_path

def isSpecies(node, nodesRank_dict):
    """This function checks if a given node is a species.
    It returns a tuple with a boolean indicating if the node is a species and the node itself."""

    if nodesRank_dict[node] == 'species':
        return True, node
    else:
        return False, node

def isGenus(node, nodesRank_dict):
    """This function checks if a given node is a genus.
    It returns a tuple with a boolean indicating if the node is a genus and the node itself."""

    if nodesRank_dict[node] == 'genus':
        return True, node
    else:
        return False, node

def isFamily(node, nodesRank_dict):
    """This function checks if a given node is a family.
    It returns a tuple with a boolean indicating if the node is a family and the node itself."""

    if nodesRank_dict[node] == 'family':
        return True, node
    else:
        return False, node
    
def isOrder(node, nodesRank_dict):
    """This function checks if a given node is an order.
    It returns a tuple with a boolean indicating if the node is an order and the node itself."""

    if nodesRank_dict[node] == 'order':
        return True, node
    else:
        return False, node

def isInfraclass(node, nodesRank_dict):
    """This function checks if a given node is an infraclass.
    It returns a tuple with a boolean indicating if the node is an infraclass and the node itself."""

    if nodesRank_dict[node] == 'infraclass':
        return True, node
    else:
        return False, node

def isSubclass(node, nodesRank_dict):
    """This function checks if a given node is a subclass.
    It returns a tuple with a boolean indicating if the node is a subclass and the node itself."""

    if nodesRank_dict[node] == 'subclass':
        return True, node
    else:
        return False, node
def isClass(node, nodesRank_dict):
    """This function checks if a given node is a class.
    It returns a tuple with a boolean indicating if the node is a class and the node itself."""

    if nodesRank_dict[node] == 'class':
        return True, node
    else:
        return False, node

def isPhylum(node, nodesRank_dict):
    """This function checks if a given node is a phylum.
    It returns a tuple with a boolean indicating if the node is a phylum and the node itself."""

    if nodesRank_dict[node] == 'phylum':
        return True, node
    else:
        return False, node
    
def isKingdom(node, nodesRank_dict):
    """This function checks if a given node is a kingdom.
    It returns a tuple with a boolean indicating if the node is a kingdom and the node itself."""

    if nodesRank_dict[node] == 'kingdom':
        return True, node
    else:
        return False, node
def isDomain(node, nodesRank_dict):
    """This function checks if a given node is a domain.
    It returns a tuple with a boolean indicating if the node is a domain and the node itself."""
    
    if nodesRank_dict[node] == 'Domain':
        return True, node
    else:
        return False, node

def getLRanks(node, nodesParent_dict, nodesRank_dict):
    """This function retrieves the lowest ranks of a given node in the taxonomy hierarchy.
    It returns a list containing the species, genus, family, order, subclass, class,
    phylum, kingdom, and domain of the node."""
    
    DFspecies = None
    DFgenus = None
    DFFamily = None
    DForder = None
    DFinfraclass = None
    DFsubclass = None
    DFclass = None
    DFphylum = None
    DFkingdom = None
    DFdomain = None

    for i in get_NodepPath(node, nodesParent_dict):
        if isSpecies(i, nodesRank_dict)[0] == True:
            DFspecies = i
        elif isGenus(i, nodesRank_dict)[0] == True:
            DFgenus = i
        elif isFamily(i, nodesRank_dict)[0] == True:
            DFFamily = i
        elif isOrder(i, nodesRank_dict)[0] == True:
            DForder = i
        elif isInfraclass(i, nodesRank_dict)[0] == True:
            DFinfraclass = i
        elif isSubclass(i, nodesRank_dict)[0] == True:
            DFsubclass = i
        elif isClass(i, nodesRank_dict)[0] == True:
            DFclass = i
        elif isPhylum(i, nodesRank_dict)[0] == True:
            DFphylum = i
        elif isKingdom(i, nodesRank_dict)[0] == True:   
            DFkingdom = i
        elif isDomain(i, nodesRank_dict)[0] == True:    
            DFdomain = i

    return [DFspecies, DFgenus, DFFamily, DForder, DFsubclass, DFinfraclass, DFclass, DFphylum, DFkingdom, DFdomain]

def getEntryRanks(blastViral, GenBank_seqs, NamesDict, NodesParent):
    """This function retrieves the ranks of entries in a GenBank sequence list.
    It returns a DataFrame containing the name, organism, tax_id, species, genus,
    family, order, subclass, class, phylum, kingdom, and domain of each entry."""
    TList = []
    for i in GenBank_seqs:
        if i.id in list(blastViral['sseqid_x']):  # this filters the sequences that are in the blast output
            if i.annotations['organism'] in NamesDict: # this filters wrong organism names in the names.dmp file
                TList.append([i.name, i.annotations['organism'], NamesDict[i.annotations['organism']]] + getLRanks(NamesDict[i.annotations['organism']], NodesParent))
    
    df = pd.DataFrame(TList, columns=['name', 'organism', 'tax_id', 'species', 'genus', 'family', 'order', 'subclass', 'class', 'phylum', 'kingdom', 'domain'])
    
    return df


def getNCBIdics(NodeFile, NamesFile):
    """    Get dictionaries from NCBI taxonomy files.
    Parameters:
    NodeFile (str): Path to the nodes.dmp file.
    NamesFile (str): Path to the names.dmp file.
    Returns:
    tuple: A tuple containing four dictionaries:
        - names_dict: Mapping of taxon names to their IDs. {'ID': 'name'}
        - IDnames_dict: Mapping of taxon IDs to their names. {'name': 'ID'}
        - nodesParent_dict: Mapping of taxon IDs to their parent IDs. {'ID': 'parent_ID'}
        - nodesRank_dict: Mapping of taxon IDs to their Linnean rank. {'ID': 'LinneanRank'}
    """
    # Read the names and nodes files into pandas DataFrames

    names_file = pd.read_csv(NamesFile, sep='\t|\t', header=None, engine='python')
    nodes_file = pd.read_csv(NodeFile, sep='\t|\t', header=None, engine='python')
    names_dict = dict(zip(names_file[2], names_file[0]))
    IDnames_dict = dict(zip(names_file[0], names_file[2]))
    nodesParent_dict = dict(zip(nodes_file[0], nodes_file[2]))
    nodesRank_dict = dict(zip(nodes_file[0], nodes_file[4]))

    return names_dict, IDnames_dict, nodesParent_dict, nodesRank_dict

def get_readCounts(kaiju_df, Lcategory):        
    """
    Get the read counts for each family in the specified category.
    
    Parameters:
    kaiju_df (DataFrame): DataFrame containing kaiju results.
    Lcategory (str): The category to filter by.
    
    Returns:
    dict: A dictionary with family names as keys and their read counts as values.
    """
    
    categoryCode = {'superkingdom': 0,
                    'phylum': 1,
                    'class': 2,
                    'order': 3,
                    'family': 4,
                    'genus': 5,
                    'species': 6}

    #taxon_list = []
    
    #for line in kaiju_df[7]:
    #    if len(line.split('; ')) > 1:
    #        taxon_list.append(line.split('; ')[categoryCode[Lcategory]])

    
    kaiju_df = kaiju_df[kaiju_df[7].isnull() == False]

    taxon_list = [line.split('; ')[categoryCode[Lcategory]] for line in kaiju_df[7]]
    
    dict_fam = {}
    
    
    for item in set(taxon_list):
        dict_fam[item] = int(0)
    
    for item in taxon_list:
        dict_fam[item] += int(1)
    
    return dict_fam