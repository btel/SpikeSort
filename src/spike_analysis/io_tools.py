#!/usr/bin/env python
#coding=utf-8
import glob, re
import os

def read_dataset(filter, dataset):
    spt = filter.read_spt(dataset)
    stim_node = "/".join(dataset.split('/')[:-1]+['stim'])
    stim = filter.read_spt(stim_node)
    return {'spt': spt['data'], 'stim': stim['data'], 'ev': []} 

def list_cells(filter, dataset):
    """List all cells which fit the pattern given in dataset. Dataset can contain
       wildcards.
       
       Example:
       
       dataset = "/Subject/sSession01/el*/cell*"
    """
    regexp = "^/(?P<subject>[a-zA-z\*]+)/s(?P<ses_id>.+)/el(?P<el_id>[0-9\*]+)/?(?P<type>[a-zA-Z]+)?(?P<cell_id>[0-9\*]+)?$"
    
    conf = filter.conf_dict
    fpath = (conf['dirname'].format(DATAPATH=os.environ['DATAPATH'])+conf['cell'])
    rec_wildcard = re.match(regexp, dataset).groupdict()
    fname = fpath.format(**rec_wildcard) 

    files = glob.glob(fname)

    rec_regexp = {"subject": "(?P<subject>[a-zA-z\*]+)",
                  "ses_id": "(?P<ses_id>.+)",
                  "el_id": "(?P<el_id>[0-9\*]+)",
                  "cell_id": "(?P<cell_id>[0-9\*]+)"}
    node_fmt = "/{subject}/s{ses_id}/el{el_id}/cell{cell_id}"
    
    f_regexp = fpath.format(**rec_regexp)
    pattern = re.compile(f_regexp)
    nodes = []
    for f in files:
        dataset_match = pattern.match(f)
        rec = rec_wildcard.copy()
        rec.update(dataset_match.groupdict())
        nodes.append(node_fmt.format(**rec))
        
    return nodes
