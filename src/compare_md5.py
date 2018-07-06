#!/usr/bin/env python3
'''
Script for comparing md5 sums of results to a reference json.
'''

import hashlib
import os
import argparse
import json


class FileWithMd5(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.basename = os.path.basename(filepath)
        self.__md5 = None

    @property
    def md5(self):
        if self.__md5 is not None:
            return self.__md5
        else:
            self.__md5 = self.calculate_md5()
            return self.__md5

    def calculate_md5(self, chunksize=4096):
        hash_md5 = hashlib.md5()
        with open(self.filepath, 'rb') as f:
            # Iter is calling f.read(chunksize) until it returns
            # the sentinel b''(empty bytes)
            for chunk in iter(lambda: f.read(chunksize), b''):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()


def get_file_with_md5(filepath):
    return FileWithMd5(filepath)


def main(args):
    with open(args.reference_json) as f:
        reference = json.load(f)
    files_to_inspect = [get_file_with_md5(file) for file in args.input_files]
    files_to_inspect_md5 = {
        file.basename: file.md5
        for file in files_to_inspect
    }
    md5_match_by_file = dict()
    match_overall = True
    try:
        for key in files_to_inspect_md5:
            match = reference[key] == files_to_inspect_md5[key]
            md5_match_by_file[key] = match
            match_overall &= match
    except KeyError:
        print('key not found')
        match_overall = False
        md5_match_by_file['match_overall'] = False
    else:
        md5_match_by_file['match_overall'] = match_overall
    with open(args.outfile, 'w') as f:
        json.dump(md5_match_by_file, fp=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', nargs='+')
    parser.add_argument('--reference_json')
    parser.add_argument('--outfile')
    args = parser.parse_args()
    main(args)
