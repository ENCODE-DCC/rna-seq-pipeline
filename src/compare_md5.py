#!/usr/bin/env python3
'''
Script for comparing md5 sums of results to a reference json.
'''

import hashlib
import os
import argparse


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
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument()
