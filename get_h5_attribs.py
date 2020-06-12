import os
import sys
import argparse
import h5py

def get_all_h5_attribs(item, levels = sys.maxsize, found_attribs = {}):
    """
    Walks the HDF5 input file/group/dataset "item" to the specified level and
    puts all the attributes it finds in the returned dict as {'/': {group:
    <{subgroups...> {dataset: {attrib: value}}}}}. Returns the dict, which will
    be empty of no attributes are found.
    """

    name = os.path.basename(item.name)
    if name == '':
        name = '/'
    #print("name:", name)

    found_attribs[name] = {}
    for attribtup in item.attrs.items():
        #print("attribtup =", attribtup)
        found_attribs[name][attribtup[0]] = attribtup[1:]

    if levels <= 1:
        return found_attribs

    levels -= 1

    if isinstance(item, h5py.Group) or isinstance(item, h5py.File):
        for key, subitem in dict(item).items():
            sub_attribs = get_all_h5_attribs(subitem, levels,
                                             found_attribs[name])

    return found_attribs

def traverse_datasets(h5_fileobj):
    '''
    Return just the list of datasets
    Modified from:
    https://stackoverflow.com/questions/51548551/reading-nested-h5-group-into-numpy-array
    '''

    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = f'{prefix}/{key}'
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    for path, _ in h5py_dataset_iterator(h5_fileobj):
        yield path

#
# Main
#
parser = argparse.ArgumentParser()
parser.add_argument(
  'hdf5_path', help='The path name of the HDF5 file'
)
parser_args = parser.parse_args()

hdf5 = h5py.File(parser_args.hdf5_path, 'r')

h5_attribs = get_all_h5_attribs(hdf5)
print('h5_attribs:', h5_attribs)

for dset in traverse_datasets(hdf5):
    print(f'<\'{dset}\', {hdf5[dset].shape}, {hdf5[dset].dtype}>,')
