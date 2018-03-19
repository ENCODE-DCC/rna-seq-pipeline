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
            indexdir='out',
            bamroot='foobar')
        aligner = align.make_aligner(star_paired)
        self.assertTrue(isinstance(aligner, align.PairedEndStarAligner))

    def test_make_aligner_star_single(self):
        star_single = argparse.Namespace(
            aligner='star',
            endedness='single',
            fastqs=['fq.fastq.gz'],
            ncpus=4,
            ramGB=8,
            indexdir='out',
            bamroot='foobar')
        aligner = align.make_aligner(star_single)
        self.assertTrue(isinstance(aligner, align.SingleEndedStarAligner))

    def create_tar_archive(self, name_type_list):
        # input : list of tuples with (name, type)
        # output: TarFile object that contains TarInfo
        # objects with specified name and type
        file_like_object = io.BytesIO(b'some content')
        archive = tarfile.open(fileobj=file_like_object, mode='w')
        for name_, type_ in name_type_list:
            tarinfo = tarfile.TarInfo(name=name_)
            tarinfo.type = type_
            archive.addfile(tarinfo)
        return archive

    def test_make_modified_TarInfo_when_one_dir_and_no_files(self):
        # setup
        archive = self.create_tar_archive([('/i/am/a/dir/', tarfile.DIRTYPE)])
        # excercise
        modified_info = align.make_modified_TarInfo(archive)
        # verify
        self.assertEqual(modified_info, [])

    def test_make_modified_TarInfo_with_one_file_and_empty_target_dir(self):
        # setup
        archive = self.create_tar_archive([('/i/am/a/file.txt',
                                            tarfile.REGTYPE)])
        # excercise
        modified_info = align.make_modified_TarInfo(archive)
        # verify
        self.assertEqual(modified_info[0].name, 'file.txt')

    def test_make_modified_TarInfo_with_one_file_one_dir_and_empty_target_dir(
            self):
        # setup
        archive = self.create_tar_archive([('/i/am/a/dir/', tarfile.DIRTYPE),
                                           ('/i/am/a/file.txt',
                                            tarfile.REGTYPE)])
        # excercise
        modified_info = align.make_modified_TarInfo(archive)
        # verify
        self.assertEqual(len(modified_info), 1)
        self.assertEqual(modified_info[0].name, 'file.txt')

    def test_make_modified_TarInfo_when_filepath_changes_to_nonempty(self):
        # setup
        archive = self.create_tar_archive([('/my/old/path/file.txt',
                                            tarfile.REGTYPE)])
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