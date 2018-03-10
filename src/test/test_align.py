#!/usr/bin/env python3
"""
Unit test module for align.py
"""

import unittest
import align
import argparse
import io
import tarfile


class TestAlignHelpers(unittest.TestCase):
    def setUp(self):
        pass

    def test_make_aligner_star_paired(self):
        star_paired = argparse.Namespace(
            aligner='star',
            endedness='paired',
            fastqs=['fq1.fastq.gz', 'fq2.fastq.gz'],
            ncpus=4,
            ramGB=8,
            indexdir='out')
        aligner = align.make_aligner(star_paired)
        self.assertTrue(isinstance(aligner, align.PairedEndStarAligner))

    def test_make_aligner_star_single(self):
        star_single = argparse.Namespace(
            aligner='star',
            endedness='single',
            fastqs=['fq.fastq.gz'],
            ncpus=4,
            ramGB=8,
            indexdir='out')
        aligner = align.make_aligner(star_single)
        self.assertTrue(isinstance(aligner, align.SingleEndedStarAligner))

    def test_make_modified_TarInfo_when_one_dir_and_no_files(self):
        '''
        fixture setup. need a hack to initialize empty archive.
        add one member that is a directory into the archive.
        '''
        file_like_object = io.BytesIO(b'some content')
        archive = tarfile.open(fileobj=file_like_object, mode='w')
        tar_directory = tarfile.TarInfo(name='/i/am/a/dir/')
        tar_directory.type = tarfile.DIRTYPE
        archive.addfile(tar_directory)
        # excercise
        modified_info = align.make_modified_TarInfo(archive)
        # verify
        self.assertEqual(modified_info, [])

    def test_make_modified_TarInfo_with_one_file_and_empty_target_dir(self):
        # setup
        file_like_object = io.BytesIO(b'some content')
        archive = tarfile.open(fileobj=file_like_object, mode='w')
        tar_file = tarfile.TarInfo(name='/i/am/a/file.txt')
        # tarfile.REGTYPE is a regular file
        tar_file.type = tarfile.REGTYPE
        archive.addfile(tar_file)
        # excercise
        modified_info = align.make_modified_TarInfo(archive)
        # verify
        self.assertEqual(modified_info[0].name, 'file.txt')

    def test_make_modified_TarInfo_with_one_file_one_dir_and_empty_target_dir(
            self):
        # setup
        file_like_object = io.BytesIO(b'some content')
        archive = tarfile.open(fileobj=file_like_object, mode='w')
        tar_file = tarfile.TarInfo(name='/i/am/a/file.txt')
        tar_file.type = tarfile.REGTYPE
        tar_directory = tarfile.TarInfo(name='/i/am/a/dir/')
        tar_directory.type = tarfile.DIRTYPE
        archive.addfile(tar_file)
        archive.addfile(tar_directory)
        # excercise
        modified_info = align.make_modified_TarInfo(archive)
        # verify
        self.assertEqual(len(modified_info), 1)
        self.assertEqual(modified_info[0].name, 'file.txt')

    def test_make_modified_TarInfo_when_filepath_changes_to_nonempty(self):
        # setup
        file_like_object = io.BytesIO(b'some content')
        archive = tarfile.open(fileobj=file_like_object, mode='w')
        tar_file = tarfile.TarInfo(name='/my/old/path/file.txt')
        tar_file.type = tarfile.REGTYPE
        archive.addfile(tar_file)
        # excercise
        modified_info = align.make_modified_TarInfo(
            archive, target_dir='/my/new/path/')
        # verify
        self.assertEqual(modified_info[0].name, '/my/new/path/file.txt')


class TestAlignerClasses(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass