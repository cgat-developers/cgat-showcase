"""rnaseq.py - Tasks for RNAseq analysis
==============================================

This module provides tasks related to RNAseq analysis.

Cufflinks/Cuffdiff
------------------

The module contains wrappers for running and parsing the output
of cufflinks and cuffdiff.

UTR analysis
------------

The functions :func:`plotGeneLevelReadExtension` and
:func:`buildUTRExtension` are part of a method to predict UTRs from
short-read data.

These methods should become part of a separate pipeline or a tool.


Requirements:

* HiddenMarkov >= 1.8.0
* cufflinks >= 2.2.1
* MASS >= 7.3.34
* RColorBrewer >= 1.0.5
* featureCounts >= 1.4.3
* samtools >= 1.1

Reference
----------

"""

import cgatcore.experiment as E
import cgatcore.csvutils as CSV
import cgat.Sra as Sra
import cgat.Fastq as Fastq

import collections
import glob
import itertools
import math
import numpy as np
import os
import pandas as pd
import re


from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.rinterface as ri

import cgat.BamTools.bamtools as BamTools
import tasks.counts as Counts
import cgatcore.database as Database
import tasks.expression as Expression
import cgat.GTF as GTF
import cgatcore.iotools as iotools
from cgatcore import pipeline as P
import tasks.mapping as mapping

from cgatcore.pipeline import cluster_runnable



class SequenceCollectionProcessor(object):
    """base class for processors of sequence collections.

    Processors of sequence collections are mappers, trimmers, etc.

    The aim of these processors is to build a sequence
    of command line statements that can be send to a single node
    on a cluster to process the input data.

    Attributes
    ==========
    compress : bool
        If True, compress temporary fastq files with gzip
    convert : bool
        If True, convert quality scores to Sanger quality scores.
    preserve_colourspace : bool
        If True, preserve colour space files. By default they are
        converted to fastq.
    tmpdir_fastq : string
        Directory with the locations of temporary :term:`fastq`
        formatted files. This directory can be used as a general
        temporary directory by a mapper.
    """

    # compress temporary fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    # set to True if you want to preserve colour space files.
    # By default, they are converted to fastq.
    preserve_colourspace = False

    # Temporary directory where temporary fastq files will
    # be located.
    tmpdir_fastq = None

    keep_sra = True

    def __init__(self, *args, **kwargs):
        pass

    def quoteFile(self, filename):
        '''return a quoted file for in-situ uncompression.

        The quoting adds uncompression for compressed files for
        programs that expect uncompressed files.

        .. note::
            This will only work if the downstream programs reads
            the file only once.

        '''
        if filename.endswith(".gz") and not self.compress:
            return "<( gunzip < %s )" % filename
        else:
            return filename

    def preprocess(self, infiles, outfile):
        '''build a preprocessing statement

        Build a command line statement that extracts/converts various
        input formats to :term:`fastq` formatted files.

        Mapping qualities are changed to solexa format.

        Returns
        -------
        statement : string
            The command line statement for pre-processing.
        files : list
            List of :term:`fastq` formatted files that will
            be created by `statement`.
        '''
        # explicitly assign qual_format to use in string formatting
        try:
            assert self.qual_format
            qual_format = self.qual_format
        except AttributeError:
            qual_format = 'illumina-1.8'

        assert len(infiles) > 0, "no input files for processing"

        tmpdir_fastq = P.get_temp_dir(shared=True)
        self.tmpdir_fastq = tmpdir_fastq
        # create temporary directory again for nodes
        statement = ["mkdir -p %s" % tmpdir_fastq]
        fastqfiles = []

        # get track by extension of outfile
        track = os.path.splitext(os.path.basename(outfile))[0]

        if self.compress:
            compress_cmd = "| gzip"
            extension = ".gz"
        else:
            compress_cmd = ""
            extension = ""

        for infile in infiles:

            if infile.endswith(".export.txt.gz"):
                # single end illumina export
                statement.append("""gunzip < %(infile)s
                | awk '$11 != "QC" || $10 ~ /(\d+):(\d+):(\d+)/ \
                {if ($1 != "")
                {readname=sprintf("%%%%s_%%%%s:%%%%s:%%%%s:%%%%s:%%%%s",
                $1, $2, $3, $4, $5, $6);}
                else {readname=sprintf("%%%%s:%%%%s:%%%%s:%%%%s:%%%%s",
                $1, $3, $4, $5, $6);}
                printf("@%%%%s\\n%%%%s\\n+\\n%%%%s\\n",readname,$9,$10);}'
                %(compress_cmd)s
                > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" % locals())
                fastqfiles.append(
                    ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))
            elif infile.endswith(".fa.gz"):
                statement.append(
                    '''gunzip < %(infile)s %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fa''' % locals())
                fastqfiles.append(("%s/%s.fa" % (tmpdir_fastq, track),))
                self.datatype = "fasta"

            elif infile.endswith(".remote"):
                files = []
                for line in iotools.open_file(infile):
                    repo, acc = line.strip().split("\t")[:2]
                    if repo == "SRA":
                        statement.append(Sra.prefetch(acc))
                        f, format, _ = Sra.peek(acc)
                        statement.append(Sra.extract(acc, tmpdir_fastq))

                        extracted_files = ["%s/%s" % (
                            tmpdir_fastq, os.path.basename(x))
                            for x in sorted(f)]
                        if not self.keep_sra:
                            statement.append(Sra.clean_cache(acc))
                        files.extend(extracted_files)

                    elif repo == "ENA":
                        filenames, dl_paths = Sra.fetch_ENA_files(acc)
                        for f in dl_paths:
                            statement.append(Sra.fetch_ENA(f, tmpdir_fastq))
                        files.extend([os.path.join(tmpdir_fastq, x) for x
                                      in filenames])

                    elif repo == "TCGA":
                        tar_name = line.strip().split("\t")[2]
                        token = glob.glob("gdc-user-token*")

                        if len(token) > 0:
                            token = token[0]
                        else:
                            token = None

                        statement.append(Sra.fetch_TCGA_fastq(acc,
                                                              tar_name,
                                                              token,
                                                              tmpdir_fastq))

                        files.append(os.path.join(
                            tmpdir_fastq, acc + "_1.fastq.gz"))
                        files.append(os.path.join(
                            tmpdir_fastq, acc + "_2.fastq.gz"))

                    else:
                        raise ValueError("Unknown repository: %s" % repo)

                if len(files) == 1 and files[0].endswith("_1.fastq.gz"):
                    new_files = [re.sub("_1.fastq.gz",
                                        ".fastq.gz",
                                        files[0])]

                elif len(files) == 2:
                    new_files = [re.sub(r"_([12]).fastq.gz",
                                        r".fastq.\1.gz",
                                        x) for x in files]

                else:
                    new_files = files

                out_base = os.path.basename(os.path.splitext(outfile)[0])
                new_files = [re.sub(r"%s/(.+).(fastq..*gz)" % tmpdir_fastq,
                                    r"%s/%s_\1.\2" % (tmpdir_fastq, out_base),
                                    nf) for nf in new_files]

                for old, new in zip(files, new_files):
                    if old != new:
                        statement.append("mv %s %s" % (old, new))

                fastqfiles.append(new_files)

            elif infile.endswith(".sra"):
                # sneak preview to determine if paired end or single end
                outdir = P.get_temp_dir()
                f, format, datatype = Sra.peek(infile)
                E.info("sra file contains the following files: %s" % f)

                # T.S need to use abi-dump for colorspace files
                if datatype == "basecalls":
                    tool = "fastq-dump"
                    self.datatype = "basecalls"
                elif datatype == "colorspace":
                    tool = "abi-dump"
                    self.datatype = "solid"

                # add extraction command to statement
                statement.append(Sra.extract(infile, tmpdir_fastq, tool))

                sra_extraction_files = ["%s/%s" % (
                    tmpdir_fastq, os.path.basename(x)) for x in sorted(f)]

                # if format is not sanger, add convert command to statement
                if 'sanger' in format and self.convert:

                    # single end fastq
                    if len(sra_extraction_files) == 1:

                        infile = sra_extraction_files[0]
                        track = os.path.splitext(os.path.basename(infile))[0]
                        if track.endswith("_1.fastq"):
                            track = track[:-8]
                        statement.append("""gunzip < %(infile)s
                        | cgat fastq2fastq
                        --method=change-format --target-format=sanger
                        --guess-format=phred64
                        --log=%(outfile)s.log
                        %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s_converted.fastq%(extension)s
                        """ % locals())

                        fastqfiles.append(
                            ("%(tmpdir_fastq)s/%(track)s_converted.fastq%(extension)s"
                             % locals(),))

                    # paired end fastqs
                    elif len(sra_extraction_files) == 2:

                        infile, infile2 = sra_extraction_files
                        track = os.path.splitext(os.path.basename(infile))[0]

                        statement.append("""gunzip < %(infile)s
                        | cgat fastq2fastq
                        --method=change-format --target-format=sanger
                        --guess-format=phred64
                        --log=%(outfile)s.log %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s_converted.1.fastq%(extension)s;
                        gunzip < %(infile2)s
                        | cgat fastq2fastq
                        --method=change-format --target-format=sanger
                        --guess-format=phred64
                        --log=%(outfile)s.log %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s_converted.2.fastq%(extension)s
                        """ % locals())

                        fastqfiles.append(
                            ("%s/%s_converted.1.fastq%s" %
                             (tmpdir_fastq, track, extension),
                             "%s/%s_converted.2.fastq%s" %
                             (tmpdir_fastq, track, extension)))
                else:
                    # Usually sra extraction files have format
                    # '1_fastq.gz' for both single and paired end files
                    # This code corrects the output to the format expected by
                    # cgat s
                    infile = sra_extraction_files[0]
                    basename = os.path.basename(infile)

                    if(len(sra_extraction_files) == 1 and
                       basename.endswith("_1.fastq.gz") and
                       self.datatype != "solid"):

                        basename = basename[:-11] + ".fastq.gz"
                        statement.append(
                            "mv %s %s/%s" % (infile, tmpdir_fastq, basename))
                        fastqfiles.append(
                            ("%s/%s" % (tmpdir_fastq, basename),))

                    # if sra extracted file is SOLID, need to keep
                    # record of qual files
                    elif self.datatype == "solid":
                        # single end SOLiD data
                        infile = P.snip(infile, "_1.fastq.gz") + \
                            "_F3.csfasta.gz"
                        quality = P.snip(
                            infile, "_F3.csfasta.gz") + "_F3_QV.qual.gz"

                        # qual file does not exist as tmpdir from
                        # SRA.extract is removed
                        # if not os.path.exists(quality):
                        #    raise ValueError("no quality file for %s" % infile)

                        fastqfiles.append((infile, quality))

                    # T.S I'm not sure if this works. Need a test case!
                    elif infile.endswith(".csfasta.F3.gz"):
                        # paired end SOLiD data
                        bn = P.snip(infile, ".csfasta.F3.gz")
                        # order is important - mirrors tophat reads followed by
                        # quals
                        f = []
                        for suffix in ("csfasta.F3", "csfasta.F5",
                                       "qual.F3", "qual.F5"):
                            fn = "%(bn)s.%(suffix)s" % locals()
                            if not os.path.exists(fn + ".gz"):
                                raise ValueError(
                                    "expected file %s.gz missing" % fn)
                            statement.append("""gunzip < %(fn)s.gz
                            %(compress_cmd)s
                            > %(tmpdir_fastq)s/%(basename)s.%(suffix)s%(extension)s
                            """ %
                                             locals())
                            f.append(
                                "%(tmpdir_fastq)s/%(basename)s.%(suffix)s%(extension)s" %
                                locals())
                        fastqfiles.append(f)

                    elif len(sra_extraction_files) == 2:
                        infile2 = sra_extraction_files[1]
                        if basename.endswith("_1.fastq.gz"):
                            basename1 = basename[:-11] + ".fastq.1.gz"
                            basename2 = basename[:-11] + ".fastq.2.gz"
                        statement.append(
                            "mv %s %s/%s; mv %s %s/%s" %
                            (infile, tmpdir_fastq, basename1,
                             infile2, tmpdir_fastq, basename2))
                        fastqfiles.append(
                            ("%s/%s" % (tmpdir_fastq, basename1),
                             "%s/%s" % (tmpdir_fastq, basename2)))

                    else:
                        fastqfiles.append(sra_extraction_files)

            elif infile.endswith(".fastq.gz"):
                format = Fastq.guessFormat(
                    iotools.open_file(infile, "r"), raises=False)
                if 'sanger' not in format and self.convert:
                    statement.append("""gunzip < %(infile)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" %
                                     locals())
                    fastqfiles.append(
                        ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))
                else:
                    E.debug("%s: assuming quality score format %s" %
                            (infile, format))
                    fastqfiles.append((infile, ))

            elif infile.endswith(".csfasta.gz"):
                # single end SOLiD data
                if self.preserve_colourspace:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"
                    if not os.path.exists(quality):
                        raise ValueError("no quality file for %s" % infile)
                    statement.append("""gunzip < %(infile)s
                    > %(tmpdir_fastq)s/%(track)s.csfasta%(extension)s""" %
                                     locals())
                    statement.append("""gunzip < %(quality)s
                    > %(tmpdir_fastq)s/%(track)s.qual%(extension)s""" %
                                     locals())
                    fastqfiles.append(("%s/%s.csfasta%s" %
                                       (tmpdir_fastq, track, extension),
                                       "%s/%s.qual%s" %
                                       (tmpdir_fastq, track, extension)))
                    self.datatype = "solid"
                else:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"

                    statement.append("""solid2fastq
                    <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)""" %
                                     locals())
                    fastqfiles.append(
                        ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))

            elif infile.endswith(".csfasta.F3.gz"):
                # paired end SOLiD data
                if self.preserve_colourspace:
                    bn = P.snip(infile, ".csfasta.F3.gz")
                    # order is important - mirrors tophat reads followed by
                    # quals
                    f = []
                    for suffix in ("csfasta.F3", "csfasta.F5",
                                   "qual.F3", "qual.F5"):
                        fn = "%(bn)s.%(suffix)s" % locals()
                        if not os.path.exists(fn + ".gz"):
                            raise ValueError(
                                "expected file %s.gz missing" % fn)
                        statement.append("""gunzip < %(fn)s.gz
                        %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s
                        """ %
                                         locals())
                        f.append(
                            "%(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s" %
                            locals())
                    fastqfiles.append(f)
                    self.datatype = "solid"
                else:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"

                    statement.append("""solid2fastq
                    <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" %
                                     locals())
                    fastqfiles.append(
                        ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))

            elif infile.endswith(".fastq.1.gz"):

                bn = P.snip(infile, ".fastq.1.gz")
                infile2 = "%s.fastq.2.gz" % bn
                if not os.path.exists(infile2):
                    raise ValueError(
                        "can not find paired ended file "
                        "'%s' for '%s'" % (infile2, infile))

                format = Fastq.guessFormat(
                    iotools.open_file(infile), raises=False)

                if 'sanger' not in format and qual_format != 'phred64':
                    statement.append("""gunzip < %(infile)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.1.fastq%(extension)s;
                    gunzip < %(infile2)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.2.fastq%(extension)s
                    """ % locals())
                    fastqfiles.append(
                        ("%s/%s.1.fastq%s" % (tmpdir_fastq, track, extension),
                         "%s/%s.2.fastq%s" % (tmpdir_fastq, track, extension)))

                elif 'sanger' not in format and qual_format == 'phred64':
                    statement.append("""gunzip < %(infile)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.1.fastq%(extension)s;
                    gunzip < %(infile2)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.2.fastq%(extension)s
                    """ % locals())
                    fastqfiles.append(
                        ("%s/%s.1.fastq%s" % (tmpdir_fastq, track, extension),
                         "%s/%s.2.fastq%s" % (tmpdir_fastq, track, extension)))
                else:
                    E.debug("%s: assuming quality score format %s" %
                            (infile, format))
                    fastqfiles.append((infile,
                                       infile2,))

            else:
                raise NotImplementedError("unknown file format %s" % infile)

        assert len(fastqfiles) > 0, "no fastq files for mapping"
        return (" ; ".join(statement) + ";", fastqfiles)


class Mapper(SequenceCollectionProcessor):
    '''Base class for short-read mappers.

    These tools map reads in :term:`fastq` or :term:`sra` formatted
    files and outut a :term:`BAM` formatted file.

    .. note::
        This class just sets the attributes below, but is not
        implementing them.

    Attributes
    ----------
    datatype : string
        Datatype of input
    strip_sequence : bool
        If True, remove sequence from BAM files in post-processing
        step.
    remove_non_unique : bool
        If True, remove non-unique matches from the BAM file. This
        processing will happen in addition to any filters applied
        in the mapper.

    Arguments
    ---------
    executable : string
        Executable to use. If unset, use default.
    strip_sequence : bool
        If True, remove sequence from BAM files in post-processing.
    remove_non_unique : bool
        If True, remove non-unique matches from the BAM file. This
        processing will happen in addition to any filters applied
        in the mapper.
    tool_options : string
        Options to be passed to the processing tool.

    '''

    datatype = "fastq"

    # strip bam files of sequenca and quality information
    strip_sequence = False

    # remove non-unique matches in a post-processing step.
    # Many aligners offer this option in the mapping stage
    # If only unique matches are required, it is better to
    # configure the aligner as removing in post-processing
    # adds to processing time.
    remove_non_unique = False

    def __init__(self,
                 executable=None,
                 strip_sequence=False,
                 remove_non_unique=False,
                 tool_options="",
                 *args, **kwargs):
        SequenceCollectionProcessor.__init__(self, *args, **kwargs)

        if executable:
            self.executable = executable
        self.strip_sequence = strip_sequence
        self.remove_non_unique = remove_non_unique

        # tool options to be passed on to the mapping tool
        self.tool_options = tool_options

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles.
        '''
        return ""

    def postprocess(self, infiles, outfile):
        '''collect output data and postprocess.'''
        return ""

    def cleanup(self, outfile):
        '''clean up.'''
        statement = '''rm -rf %s;''' % (self.tmpdir_fastq)

        return statement

    def build(self, infiles, outfile):
        '''run mapper

        This method combines the output of the :meth:`preprocess`,
        :meth:`mapper`, :meth:`postprocess` and :meth:`clean` sections
        into a single statement.

        Arguments
        ---------
        infiles : list
             List of input filenames
        outfile : string
             Output filename

        Returns
        -------
        statement : string
             A command line statement. The statement can be a series
             of commands separated by ``;`` and/or can be unix pipes.

        '''

        cmd_preprocess, mapfiles = self.preprocess(infiles, outfile)
        cmd_mapper = self.mapper(mapfiles, outfile)
        cmd_postprocess = self.postprocess(infiles, outfile)
        cmd_clean = self.cleanup(outfile)

        assert cmd_preprocess.strip().endswith(";"),\
            "missing ';' at end of command %s" % cmd_preprocess.strip()
        assert cmd_mapper.strip().endswith(";"),\
            "missing ';' at end of command %s" % cmd_mapper.strip()
        if cmd_postprocess:
            assert cmd_postprocess.strip().endswith(";"),\
                "missing ';' at end of command %s" % cmd_postprocess.strip()
        if cmd_clean:
            assert cmd_clean.strip().endswith(";"),\
                "missing ';' at end of command %s" % cmd_clean.strip()

        statement = " ".join((cmd_preprocess,
                              cmd_mapper,
                              cmd_postprocess,
                              cmd_clean))

        return statement


class Kallisto(Mapper):

    '''run Kallisto to quantify transcript abundance from fastq files
    - set pseudobam to True to output a pseudobam along with the quantification'''

    def __init__(self, pseudobam=False, readable_suffix=False, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

        self.pseudobam = pseudobam
        self.readable_suffix = readable_suffix

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles'''

        tmpdir = os.path.join(self.tmpdir_fastq + "kallisto")

        logfile = outfile + ".log"

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        if nfiles == 1:

            infiles = (" --fragment-length=%(kallisto_fragment_length)s" +
                       " --sd=%(kallisto_fragment_sd)s" +
                       " --single %s" % " ".join([x[0] for x in infiles]))

        elif nfiles == 2:
            infiles = " ".join([" ".join(x) for x in infiles])

        else:
            raise ValueError("incorrect number of input files")

        outdir = os.path.dirname(outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        statement = '''
        kallisto quant %%(kallisto_options)s
        -t %%(job_threads)s
        --bootstrap-samples=%%(kallisto_bootstrap)s
        -i %%(index)s -o %(tmpdir)s %(infiles)s''' % locals()

        if self.pseudobam:
            statement += '''
            --pseudobam | samtools view -b -
            > %(outfile)s.bam 2> %(logfile)s;''' % locals()
        else:
            statement += ''' > %(logfile)s &> %(logfile)s ;''' % locals()

        self.tmpdir = tmpdir

        return statement

    def postprocess(self, infiles, outfile):
        '''move outfiles from tmpdir to final location'''

        tmpdir = self.tmpdir

        statement = ('''
        mv -f %(tmpdir)s/abundance.h5 %(outfile)s;
        ''' % locals())

        if self.readable_suffix:

            outfile_readable = outfile + self.readable_suffix
            statement += ('''
            kallisto h5dump -o %(tmpdir)s %(outfile)s;
            mv %(tmpdir)s/abundance.tsv %(outfile_readable)s;
            rm -rf %(tmpdir)s/bs_abundance_*.tsv;''' % locals())

        return statement

    def cleanup(self, outfile):
        '''clean up.'''
        statement = '''rm -rf %s; rm -rf %s;''' % (
            self.tmpdir_fastq, self.tmpdir)

        return statement
        
        
def normaliseCounts(counts_inf, outfile, method):
    ''' normalise a counts table'''

    counts = Counts.Counts(counts_inf)
    counts.normalise(method=method)
    counts.table.index.name = "id"

    if outfile.endswith(".gz"):
        counts.table.to_csv(outfile, sep="\t", compression="gzip")
    else:
        counts.table.to_csv(outfile, sep="\t")


def parse_table(sample, outfile_raw, outfile, columnname):
    '''
    parse the output of featurecounts or alignment free qauntifiers
    and extract number of reads for downstream quantification
    '''

    column_ix = findColumnPosition(outfile_raw, columnname)
    sample = sample

    if outfile_raw.endswith("gz"):
        grep = "zgrep"
    else:
        grep = "grep"

    statement = '''
    echo -e "id\\t%(sample)s" | gzip > %(outfile)s;
    %(grep)s -v '^#' %(outfile_raw)s |
    cut -f1,%(column_ix)s | awk 'NR>1' | gzip >> %(outfile)s;
    ''' % locals()
    P.run(statement)


def estimateSleuthMemory(bootstraps, samples, transcripts):
    ''' The memory usage of Sleuth is dependent upon the number of
    samples, transcripts and bootsraps.

    A rough estimate is:
    24 bytes * bootstraps * samples * transcripts
    (https://groups.google.com/forum/#!topic/kallisto-sleuth-users/mp064J-DRfI)

    TS: I've found this to be a serious underestimate so we use a
    more conservative estimate here (48 bytes * ... ) with a default of 2G
    '''

    estimate = (48 * bootstraps * samples * transcripts)

    job_memory = "%fG" % (max(2.0, (estimate / 1073741824.0)))

    return job_memory


def findColumnPosition(infile, column):
    ''' find the position in the header of the specified column
    The returned value is one-based (e.g for bash cut)'''
    with iotools.open_file(infile, "r") as inf:
        while True:
            header = inf.readline()
            if not header.startswith("#"):
                head = header.strip().split("\t")
                j = 0
                for h in head:
                    if column in h:
                        column_ix = j
                    j += 1
                # column_ix = header.strip().split("\t").index(column)
                break
        if column_ix:
            return column_ix + 1
        else:
            raise ValueError("could not find %s in file header" % column)


class Quantifier(object):
    ''' base class for transcript and gene-level quantification from a
    BAM or fastq

    Note: All quantifier classes are designed to perform both
    transcript-level and gene-level analyses in a single run. The
    runGene method
    '''

    def __init__(self, infile, transcript_outfile, gene_outfile, annotations,
                 job_threads=None, job_memory=None, strand=None,
                 options=None, bootstrap=None,
                 fragment_length=None, fragment_sd=None,
                 transcript2geneMap=None, libtype=None, kmer=None,
                 biascorrect=None):
        '''
        Attributes
        ----------
        infile: string
           Input  filename
        transcript_outfile: string
           Outfile of transcript quantifications in :term: `gz.raw` format
        gene_outfile: string
           Outfile of gene quantifications in :term: `gz.raw` format
        job_threads: string
           Number of threads per job
        strand: int
           For FeatureCounts the strand is specified as either 0, 1, 2
        options: string
           Options specified as a string
        annotations: string
           Filename with gene set in :term:`gtf` format.
        bootstrap: int
           Number of boostrap values for alignment free quantifiers
        job_memory: str
           Amount of memory available for job
        frangment_length: int
           Must specify the expected fragment length for single-end reads
           This is specified in pipeline_ini.
           :term:`PARAMS` - fragment_length option.
        frangment_sd: int
           Must specify the expected fragment length sd for single-end reads
           This is specified in pipeline_ini.
           :term:`PARAMS` - fragment_sd option.
        libtype: string
           This is specified in pipeline_ini
           :term:`PARAMS` - library type option.
        kmer: int
           This is specified in the pipeline.yml
           :term:`PARAMS` - kmer size for aligment free quant.
        '''

        self.infile = infile
        self.transcript_outfile = transcript_outfile
        self.gene_outfile = gene_outfile
        self.job_threads = job_threads
        self.strand = strand
        self.options = options
        self.annotations = annotations
        self.bootstrap = bootstrap
        self.job_memory = job_memory
        self.fragment_length = fragment_length
        self.fragment_sd = fragment_sd
        self.t2gMap = transcript2geneMap
        self.libtype = libtype
        self.biascorrect = biascorrect

        # TS: assume sample name is directory for outfile which is
        # should be for pipeline_rnaseqdiffexpression. This would be
        # better handled in pipeline though
        self.sample = os.path.basename(os.path.dirname(self.gene_outfile))

    def run_transcript(self):
        ''' generate transcript-level quantification estimates'''

    def run_gene(self):
        ''' generate gene-level quantification estimates'''

    def run_all(self):
        ''' '''
        self.run_transcript()
        self.run_gene()


class AF_Quantifier(Quantifier):
    ''' Parent class for all alignment-free quantification methods'''

    def run_gene(self):
        ''' Aggregate transcript counts to generate gene-level counts
        using a map of transript_id to gene_id '''

        transcript_df = pd.read_table(self.transcript_outfile,
                                      sep="\t", index_col=0)
        transcript2gene_df = pd.read_table(self.t2gMap, sep="\t", index_col=0)
        transcript_df = pd.merge(transcript_df, transcript2gene_df,
                                 left_index=True, right_index=True,
                                 how="inner")
        gene_df = pd.DataFrame(transcript_df.groupby(
             'gene_id')[self.sample].sum())
        gene_df.index.name = 'id'

        gene_df.to_csv(
            self.gene_outfile, compression="gzip", sep="\t")


class KallistoQuantifier(AF_Quantifier):
    ''' quantifier class to run kallisto'''

    def run_transcript(self):
        ''' '''
        fastqfile = self.infile
        index = self.annotations
        job_threads = self.job_threads
        job_memory = self.job_memory
        kallisto_options = self.options
        kallisto_bootstrap = self.bootstrap
        kallisto_fragment_length = self.fragment_length
        kallisto_fragment_sd = self.fragment_sd
        outfile = os.path.join(
            os.path.dirname(self.transcript_outfile), "abundance.h5")
        sample = self.sample

        # kallisto output is in binary (".h5") format
        # Supplying a "readable_suffix" to the mapping.Kallisto
        # ensures an additional human readable file is also generated
        readable_suffix = ".tsv"
        m = Kallisto(readable_suffix=readable_suffix)

        statement = m.build((fastqfile), outfile)

        P.run(statement)

        outfile_readable = outfile + readable_suffix

        # parse the output to extract the counts
        parse_table(self.sample, outfile_readable,
                    self.transcript_outfile, 'est_counts')


@cluster_runnable
def makeExpressionSummaryPlots(counts_inf, design_inf, logfile):
    ''' use the plotting methods for Counts object to make summary plots'''

    with iotools.open_file(logfile, "w") as log:

        plot_prefix = P.snip(logfile, ".log")
        log.write("1")

        # need to manually read in data as index column is not the first column
        in_table = pd.read_table(counts_inf, sep="\t", index_col=0)
        in_table = in_table.dropna(axis=0)
        counts = Counts.Counts(in_table)
        counts.table.columns = [x.replace(".", "-") for x in counts.table.columns]
        log.write("2")
        design = Expression.ExperimentalDesign(design_inf)
        log.write("3")
        # make certain counts table only include samples in design
        counts.restrict(design)
        log.write("4")
        cor_scatter_outfile = plot_prefix + "_pairwise_correlations_scatter.png"
        cor_heatmap_outfile = plot_prefix + "_pairwise_correlations_heatmap.png"
        pca_var_outfile = plot_prefix + "_pca_variance.png"
        pca1_outfile = plot_prefix + "_pc1_pc2.png"
        pca2_outfile = plot_prefix + "_pc3_pc4.png"
        heatmap_outfile = plot_prefix + "_heatmap.png"

        # use log expression so that the PCA is not overly biased
        # towards the variance in the most highly expressed genes
        counts_log10 = counts.log(base=10, pseudocount=0.1, inplace=False)
        log.write("2")
        log.write("plot correlations scatter: %s\n" % cor_scatter_outfile)
        log.write("plot correlations heatmap: %s\n" % cor_heatmap_outfile)
        # counts_log10.plotPairwise(
        # cor_scatter_outfile, cor_heatmap_outfile, subset=2000)

        # for the heatmap, and pca we want the top expressed genes (top 25%).
        counts_log10.removeObservationsPerc(percentile_rowsums=75)

        log.write("plot pc1,pc2: %s\n" % pca1_outfile)
        counts_log10.plotPCA(design,
                             pca_var_outfile, pca1_outfile,
                             x_axis="PC1", y_axis="PC2",
                             colour="group", shape="group")

        log.write("plot pc3,pc4: %s\n" % pca2_outfile)
        counts_log10.plotPCA(design,
                             pca_var_outfile, pca2_outfile,
                             x_axis="PC3", y_axis="PC4",
                             colour="group", shape="group")

        # Z-score normalise the expression for the heatmap visualisation
        # counts_log10.zNormalise(inplace=True)

        # log.write("plot heatmap: %s\n" % heatmap_outfile)
        # counts_log10.heatmap(heatmap_outfile, zscore=True)


def getAlignmentFreeNormExp(transcript_infiles, basename, column,
                            transcripts_outf, genes_outf, t2gMap):
    ''' Extract the normalised expression from the transcript-level
    quantification, merge across multiple samples and output
    transcript-level and gene-level tables'''

    for n, infile in enumerate(transcript_infiles):
        # replace filename to use the full results table
        dirname = os.path.dirname(infile)
        sample = os.path.basename(dirname)
        infile = os.path.join(dirname, basename)

        tmp_df = pd.read_table(infile, sep="\t", index_col=0)
        tmp_df.drop([x for x in tmp_df.columns if x != column],
                    axis=1, inplace=True)
        tmp_df.columns = [sample]
        tmp_df.index.name = "id"

        if n == 0:
            transcript_df = tmp_df
        else:
            transcript_df = transcript_df.merge(
                tmp_df, how="outer", left_index=True, right_index=True)

    transcript_df.to_csv(transcripts_outf, sep="\t", compression="gzip")

    transcript2gene_df = pd.read_table(t2gMap, sep="\t", index_col=0)
    transcript_df = pd.merge(transcript_df, transcript2gene_df,
                             left_index=True, right_index=True,
                             how="inner")
    gene_df = pd.DataFrame(transcript_df.groupby('gene_id').sum())
    gene_df.index.name = 'id'

    gene_df.to_csv(
        genes_outf, compression="gzip", sep="\t")

'''
# ########## old code ################
'''


def filterAndMergeGTF(infile, outfile, remove_genes, merge=False):
    '''remove genes from GTF file.

    Genes that match gene identifiers the dictionary `remove_genes`
    are removed. The dictionary maps lists of labels to genes. The
    labels can for example correspond to filters that have been
    applied and that require this particular gene to be removed, for
    example::

        remove_genes = {'ENSG00001' : ('is_repeat', 'is_known'),
                        'ENSG00003' : ('is_repeat')}

    A summary file :file:`<outfile>.summary` contains the number of
    transcripts that failed various filters.

    The :file:`<outfile>.removed.tsv.gz` contains the filters that a
    transcript failed.

    Arguments
    ---------
    infile : string
        Input filename in :term:`gtf` format.
    outfile : string
        Output filename in :term:`gtf` format.
    remove_genes : dict
        Dictionary mapping gene identifiers to names of filters.
    merge : bool
        If True, the resultant transcript models are merged by overlap.

    '''

    counter = E.Counter()

    # write summary table
    outf = iotools.open_file(outfile + ".removed.tsv.gz", "w")
    outf.write("gene_id\tnoverlap\tsection\n")
    for gene_id, r in remove_genes.items():
        for s in r:
            counter[s] += 1
        outf.write("%s\t%i\t%s\n" % (gene_id,
                                     len(r),
                                     ",".join(r)))
    outf.close()

    # filter gtf file
    tmpfile = P.get_temp_file(".")
    inf = GTF.iterator(iotools.open_file(infile))

    genes_input, genes_output = set(), set()

    for gtf in inf:
        genes_input.add(gtf.gene_id)
        if gtf.gene_id in remove_genes:
            continue
        genes_output.add(gtf.gene_id)
        tmpfile.write("%s\n" % str(gtf))

    tmpfile.close()
    tmpfilename = tmpfile.name

    outf = iotools.open_file(outfile + ".summary.tsv.gz", "w")
    outf.write("category\ttranscripts\n")
    for x, y in counter.items():
        outf.write("%s\t%i\n" % (x, y))
    outf.write("input\t%i\n" % len(genes_input))
    outf.write("output\t%i\n" % len(genes_output))
    outf.write("removed\t%i\n" % (len(genes_input) - len(genes_output)))

    outf.close()

    # close-by exons need to be merged, otherwise
    # cuffdiff fails for those on "." strand

    if merge:
        statement = '''
        sort -t$'\\t' -k1,1 -k4,4n
        < %(tmpfilename)s
        | cgat gtf2gtf
            --method=unset-genes --pattern-identifier="NONC%%06i"
            --log=%(outfile)s.log
        | cgat gtf2gtf
            --method=merge-genes
            --log=%(outfile)s.log
        | cgat gtf2gtf
            --method=merge-exons
            --merge-exons-distance=5
            --log=%(outfile)s.log
        | cgat gtf2gtf
            --method=renumber-genes
            --pattern-identifier="NONC%%06i"
            --log=%(outfile)s.log
        | cgat gtf2gtf
            --method=renumber-transcripts
            --pattern-identifier="NONC%%06i"
            --log=%(outfile)s.log
        | sort -t$'\\t' -k1,1 -k9,9 -k4,4n
        | gzip > %(outfile)s
        '''
    else:
        statement = '''
        sort -t$'\\t' -k1,1 -k4,4n
        < %(tmpfilename)s
        | gzip > %(outfile)s
        '''

    P.run(statement)

    os.unlink(tmpfilename)


def runCufflinks(gtffile, bamfile, outfile, job_threads=1):
    '''run cufflinks to estimate expression levels.

    See cufflinks manuals for full explanation of infiles/outfiles/options
    http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html

    Arguments
    ---------
    gtffile : string
        Filename of geneset in :term:`gtf` format.

    bamfile : string
        Filename of reads in :term:`bam` format.

    genome_dir : string
        :term:`PARAMS` - genome directory containing fasta file. This is
        specified in pipeline_ini

    cufflinks_library_type : string
        :term:`PARAMS` - cufflinks library type option. This is
        specified in pipeline_ini

    cufflinks_options : string
        :term:`PARAMS` - cufflinks options (see manual). These are
        specified in pipeline_ini

    outfile : string
        defines naming of 3 output files for each input file

        1.outfile.gtf.gz:  transcripts.gtf file in :term:`gtf` format
        produced by cufflinks (see manual). Contains the assembled gene
        isoforms.
        This is the file used for the downstream file analysis

        2.outfile.fpkm_tracking.gz: renamed outfile.isoforms.fpkm_tracking file
        from cufflinks - contains estimated isoform-level
        expression values in "FPKM Tracking Format".

        3.outfile.genes_tracking.gz: renamed outfile.genes.fpkm_tracking.gz
        from cufflinks - contains estimated gene-level
        expression values in "FPKM Tracking Format".

    job_threads : int
        Number of threads to use
    '''

    track = os.path.basename(P.snip(gtffile, ".gtf.gz"))

    tmpdir = P.get_temp_dir()

    gtffile = os.path.abspath(gtffile)
    bamfile = os.path.abspath(bamfile)
    outfile = os.path.abspath(outfile)

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    # increase max-bundle-length to 4.5Mb due to Galnt-2 in mm9 with a
    # 4.3Mb intron.

    # AH: removed log messages about BAM record error
    # These cause logfiles to grow several Gigs and are
    # frequent for BAM files not created by tophat.

    # Error is:
    # BAM record error: found spliced alignment without XS attribute
    statement = '''mkdir %(tmpdir)s;
    cd %(tmpdir)s;
    cufflinks --label %(track)s
              --GTF <(gunzip < %(gtffile)s)
              --num-threads %(job_threads)i
              --frag-bias-correct %(genome_dir)s/%(genome)s.fa
              --library-type %(cufflinks_library_type)s
              %(cufflinks_options)s
              %(bamfile)s
    | grep -v 'BAM record error'
    >& %(outfile)s;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s.gtf.gz;
    gzip < isoforms.fpkm_tracking > %(outfile)s.fpkm_tracking.gz;
    gzip < genes.fpkm_tracking > %(outfile)s.genes_tracking.gz;
    rm -rf %(tmpdir)s
    '''

    P.run(statement)


def loadCufflinks(infile, outfile):
    '''load cufflinks expression levels into database

        Takes cufflinks output and loads into database for later report
        building
        For each input file it generates two tables in a sqlite database:

        1. outfile_fpkm: contains information from infile.fpkm_tracking.gz
        2. outfile_genefpkm : contains information from
        infile.genes_tracking.gz

    Arguments
    ---------
    infile : string
        Cufflinks output. This is used to find
        auxiliary files: specifically infile.genes_tracking.gz and
        infile.fpkm_tracking.gz
    outfile : string
        Output filename used to create logging information in `.load` files.
        Also used to create "_fpkm" and "_genefpkm" tables in database.
    '''

    track = P.snip(outfile, ".load")
    P.load(infile + ".genes_tracking.gz",
           outfile=track + "_genefpkm.load",
           options="--add-index=gene_id "
           "--ignore-column=tracking_id "
           "--ignore-column=class_code "
           "--ignore-column=nearest_ref_id")

    track = P.snip(outfile, ".load")
    P.load(infile + ".fpkm_tracking.gz",
           outfile=track + "_fpkm.load",
           options="--add-index=tracking_id "
           "--ignore-column=nearest_ref_id "
           "--rename-column=tracking_id:transcript_id")

    iotools.touch_file(outfile)


def quantifyWithStringTie(gtffile, bamfile, outdir):
    '''Run string tie in quantitation mode

    Arguments
    ---------
    gtffile: string
        Name of gtf/gtf.gz file to quantify against
    bamfile: string
        Name of BAM file with alignments to quantify with
    outdir: string
        Directory to place output files in

    The output files generated are:
    1. e_data.ctab: Exon-level expression measurements
    2. i_data.ctab: intron/junction level expression measurements
    3. t_data.ctab: transcript level expression measurements
    4. e2t.ctab: exon 2 transcript lookup table
    5. i2t.ctab: intron to transcript lookup table

    Parameters
    ----------

    quant_threads: int
        Number of threads to use
    quant_options: string
        Commandline arguements to StringTie
    quant_memory: string
        Memory to request for job'''

    if gtffile.endswith(".gz"):
        gtffile = "<( zcat %s)" % gtffile

    job_threads = PARAMS["stringtie_quant_threads"]
    job_memory = PARAMS["stringtie_quant_memory"]

    statement = '''stringtie -G %(gtffile)s
                                %(bamfile)s
                             -eb %(outdir)s
                             %(stringtie_quant_options)s
                             > /dev/null
                             2> %(outdir)s.log'''

    if bamfile.endswith(".remote"):
        token = glob.glob("gdc-user-token*")
        tmpfilename = P.get_temp_filename()
        if len(token) > 0:
            token = token[0]
        else:
            token = None

        s, bamfile = Sra.process_remote_BAM(
            bamfile, token, tmpfilename,
            filter_bed=os.path.join(
                PARAMS["annotations_dir"],
                PARAMS["annotations_interface_contigs_bed"]))

        bamfile = " ".join(bamfile)
        statement = " ; ".join(
            ["mkdir %(tmpfilename)s",
             s,
             statement,
             "rm -r %(tmpfilename)s"])

    P.run(statement)


def merge_and_loadStringTie(infiles, track_regex, outfile):
    '''Load stringtie quantitation from multiple tracks into a set
    of database tables.

    Arguements
    ----------
    infiles: list of list of string
        infiles contains the output tables from a stringtie -b run for
        several tracks. Each item in the the list is the output files
        from a single track. There are five output files (see
        :function:rnaseq.quantifyWithStringTie). It is assumed
        the files within each list are in the same order in a
        directory named after the track. e.g.
        [["track1/i_data.ctab", ...],["track2/i_data.ctab",...],...]
    track_regex: string
        regular expression to capture track name out of filenames
    outfile: string
        output file, should end in .load, is used for table prefix

    Adds the following tables to the database:
        PREFIX_transcript_data
        PREFIX_exon_data
        PREFIX_intron_data
        PREFIX_exon2transcript
        PREFIX_intron2transcript

    The first three will have a track column indicating which track
    they come from.  The last two are assumed to be identical for each
    track (this is tested and confirmed) '''

    infiles = zip(*infiles)

    table_suffixes = {"i_data.ctab": "intron_data",
                      "e_data.ctab": "exon_data",
                      "t_data.ctab": "transcript_data",
                      "e2t.ctab": "exon2transcript",
                      "i2t.ctab": "intron2transcript"}

    table_indexes = {"i_data.ctab": "i_id",
                     "e_data.ctab": "e_id",
                     "t_data.ctab": "t_id,t_name,gene_id",
                     "e2t.ctab": "e_id,t_id",
                     "i2t.ctab": "i_id,t_id"}

    table_prefix = P.snip(outfile, ".load")

    for infile in infiles:
        which_files = set([os.path.basename(f) for f in infile])
        assert len(which_files) == 1, "Input file lists not in same order"
        table_suffix = table_suffixes[list(which_files)[0]]

        tablename = os.path.basename(table_prefix + "_" + table_suffix)
        indexs = table_indexes[list(which_files)[0]]

        if "2" in table_suffix:
            P.load(infile[0], outfile,
                   options="-i %s" % indexs,
                   tablename=tablename)
            continue

        P.concatenate_and_load(infile, outfile,
                               regex_filename=track_regex,
                               tablename=tablename,
                               options="--quick -i %s" % indexs,
                               job_memory="12G")


def mergeCufflinksFPKM(infiles, outfile, genesets,
                       tracking="genes_tracking",
                       identifier="gene_id"):

    '''build aggregate table with cufflinks FPKM values.


    Arguments
    ---------
    infiles : list
        Filenames of cufflinks results
    outfile : string
        Output filename in :term:`tsv` format.
    genesets : string
        Genesets that have been used. This is used
        to derive the prefix for extracting the correct
        column from the cufflinks output.
    tracking : string
        Select file type to merge. Valid values are `genes_tracking`
        to merge per-gene estimates and `fpkm_tracking` to merge
        per-transcript estimates.
    identifier : string
        Identifier to use for genes/transcripts. Replaces the
        default `tracking_id` from cufflinks.
    '''

    for x in genesets:
        if str(x) in os.path.basename(outfile):
            prefix = str(x)

    headers = ",".join(
        [re.match("fpkm.dir/.*_(.*).cufflinks", x).groups()[0]
         for x in infiles])

    statement = '''
    cgat combine_tables
        --log=%(outfile)s.log
        --columns=1
        --skip-titles
        --header-names=%(headers)s
        --take=FPKM fpkm.dir/%(prefix)s_*.%(tracking)s.gz
    | perl -p -e "s/tracking_id/%(identifier)s/"
    | (read h; echo \"$h\"; sort ) "
    | gzip
    > %(outfile)s
    '''
    P.run(statement)


def runFeatureCounts(annotations_file,
                     bamfile,
                     outfile,
                     job_threads=4,
                     strand=0,
                     options=""):
    '''run FeatureCounts to collect read counts.

    If `bamfile` is paired, paired-end counting is enabled and the bam
    file automatically sorted.

    Arguments
    ---------
    annotations_file : string
        Filename with gene set in :term:`gtf` format.
    bamfile : string
        Filename with short reads in :term:`bam` format.
    outfile : string
        Output filename in :term:`tsv` format.
    job_threads : int
        Number of threads to use.
    strand : int
        Strand option in FeatureCounts.
    options : string
        Options to pass on to FeatureCounts.

    '''

    # featureCounts cannot handle gzipped in or out files
    outfile = P.snip(outfile, ".gz")
    tmpdir = P.get_temp_dir()
    annotations_tmp = os.path.join(tmpdir,
                                   'geneset.gtf')
    bam_tmp = os.path.join(tmpdir,
                           os.path.basename(bamfile))

    # -p -B specifies count fragments rather than reads, and both
    # reads must map to the feature
    # for legacy reasons look at feature_counts_paired
    if BamTools.is_paired(bamfile):
        # select paired end mode, additional options
        paired_options = "-p -B"
        # sort by read name
        paired_processing = \
            """samtools
            sort -@ %(job_threads)i -n -o %(bam_tmp)s %(bamfile)s;
            """ % locals()
        bamfile = bam_tmp
    else:
        paired_options = ""
        paired_processing = ""

    # AH: what is the -b option doing?
    statement = '''mkdir %(tmpdir)s;
                   zcat %(annotations_file)s > %(annotations_tmp)s;
                   %(paired_processing)s
                   featureCounts %(options)s
                                 -T %(job_threads)i
                                 -s %(strand)s
                                 -a %(annotations_tmp)s
                                 %(paired_options)s
                                 -o %(outfile)s
                                 %(bamfile)s
                    >& %(outfile)s.log;
                    gzip -f %(outfile)s;
                    rm -rf %(tmpdir)s
    '''

    P.run(statement)


def buildExpressionStats(
        dbhandle,
        outfile,
        tablenames,
        outdir,
        regex_table="(?P<design>[^_]+)_"
        "(?P<geneset>[^_]+)_"
        "(?P<counting_method>[^_]+)_"
        "(?P<method>[^_]+)_"
        "(?P<level>[^_]+)_diff"):
    """compile expression summary statistics from database.

    This method outputs a table with the number of genes tested,
    failed, differentially expressed, etc. for a series of DE tests.

    Arguments
    ---------
    dbhandle : object
        Database handle.
    tables : list
        List of tables to process.
    outfile : string
        Output filename in :term:`tsv` format.
    outdir : string
        Output directory for diagnostic plots.
    regex : string
        Regular expression to extract experimental information
        from table name.
    """

    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = iotools.open_file(outfile, "w")
    outf.write("\t".join(
        ("design",
         "geneset",
         "level",
         "counting_method",
         "treatment_name",
         "control_name",
         "tested",
         "\t".join(["status_%s" % x for x in keys_status]),
         "significant",
         "twofold")) + "\n")

    for tablename in tablenames:
        r = re.search(regex_table, tablename)
        if r is None:
            raise ValueError(
                "can't match tablename '%s' to regex" % tablename)
        geneset = r.group("geneset")
        design = r.group("design")
        level = r.group("level")
        counting_method = r.group("counting_method")
        geneset = r.group("geneset")

        def toDict(vals, l=2):
            return collections.defaultdict(
                int,
                [(tuple(x[:l]), x[l]) for x in vals])

        tested = toDict(Database.executewait(
            dbhandle,
            "SELECT treatment_name, control_name, "
            "COUNT(*) FROM %(tablename)s "
            "GROUP BY treatment_name,control_name" % locals()
            ).fetchall())
        status = toDict(Database.executewait(
            dbhandle,
            "SELECT treatment_name, control_name, status, "
            "COUNT(*) FROM %(tablename)s "
            "GROUP BY treatment_name,control_name,status"
            % locals()).fetchall(), 3)
        signif = toDict(Database.executewait(
            dbhandle,
            "SELECT treatment_name, control_name, "
            "COUNT(*) FROM %(tablename)s "
            "WHERE significant "
            "GROUP BY treatment_name,control_name" % locals()
            ).fetchall())

        fold2 = toDict(Database.executewait(
            dbhandle,
            "SELECT treatment_name, control_name, "
            "COUNT(*) FROM %(tablename)s "
            "WHERE (l2fold >= 1 or l2fold <= -1) AND significant "
            "GROUP BY treatment_name,control_name,significant"
            % locals()).fetchall())

        for treatment_name, control_name in tested.keys():
            outf.write("\t".join(map(str, (
                design,
                geneset,
                level,
                counting_method,
                treatment_name,
                control_name,
                tested[(treatment_name, control_name)],
                "\t".join(
                    [str(status[(treatment_name, control_name, x)])
                     for x in keys_status]),
                signif[(treatment_name, control_name)],
                fold2[(treatment_name, control_name)]))) + "\n")

        # plot length versus P-Value
        data = Database.executewait(
            dbhandle,
            "SELECT i.sum, pvalue "
            "FROM %(tablename)s, "
            "%(geneset)s_geneinfo as i "
            "WHERE i.gene_id = test_id AND "
            "significant" % locals()).fetchall()

        # require at least 10 datapoints - otherwise smooth scatter fails
        if len(data) > 10:
            data = zip(*data)

            pngfile = ("%(outdir)s/%(design)s_%(geneset)s_%(level)s"
                       "_pvalue_vs_length.png") % locals()
            R.png(pngfile)
            R.smoothScatter(R.log10(ro.FloatVector(data[0])),
                            R.log10(ro.FloatVector(data[1])),
                            xlab='log10( length )',
                            ylab='log10( pvalue )',
                            log="x", pch=20, cex=.1)

            R['dev.off']()

    outf.close()


def loadCuffdiff(dbhandle, infile, outfile, min_fpkm=1.0):
    '''load results from cuffdiff analysis to database

    This functions parses and loads the results of a cuffdiff differential
    expression analysis.
    Parsing is performed by the parseCuffdiff function.

    Multiple tables will be created as cuffdiff outputs information
    on gene, isoform, tss, etc. levels.

    The method converts from ln(fold change) to log2 fold change.

    Pairwise comparisons in which one gene is not expressed (fpkm <
    `min_fpkm`) are set to status 'NOCALL'. These transcripts might
    nevertheless be significant.

    Arguments
    ---------
    dbhandle : object
        Database handle.
    infile : string
        Input filename, output from cuffdiff
    outfile : string
        Output filename in :term:`tsv` format.
    min_fpkm : float
        Minimum fpkm. Genes with an fpkm lower than this will
        be set to status `NOCALL`.

    '''

    prefix = P.to_table(outfile)
    indir = infile + ".dir"

    if not os.path.exists(indir):
        iotools.touch_file(outfile)
        return

    # E.info( "building cummeRbund database" )
    # R('''library(cummeRbund)''')
    # cuff = R('''readCufflinks(dir = %(indir)s, dbfile=%(indir)s/csvdb)''' )
    # to be continued...

    tmpname = P.get_temp_filename(shared=True)

    # ignore promoters and splicing - no fold change column, but  sqrt(JS)
    for fn, level in (("cds_exp.diff.gz", "cds"),
                      ("gene_exp.diff.gz", "gene"),
                      ("isoform_exp.diff.gz", "isoform"),
                      # ("promoters.diff.gz", "promotor"),
                      # ("splicing.diff.gz", "splice"),
                      ("tss_group_exp.diff.gz", "tss")):

        tablename = prefix + "_" + level + "_diff"

        infile = os.path.join(indir, fn)

        results = parseCuffdiff(infile, min_fpkm=min_fpkm)
        Expression.writeExpressionResults(tmpname, results)
        P.load(tmpname, outfile,
               tablename=tablename,
               options="--allow-empty-file "
               "--add-index=treatment_name "
               "--add-index=control_name "
               "--add-index=test_id")

    for fn, level in (("cds.fpkm_tracking.gz", "cds"),
                      ("genes.fpkm_tracking.gz", "gene"),
                      ("isoforms.fpkm_tracking.gz", "isoform"),
                      ("tss_groups.fpkm_tracking.gz", "tss")):

        tablename = prefix + "_" + level + "_levels"
        infile = os.path.join(indir, fn)

        P.load(infile, outfile,
               tablename=tablename,
               options="--allow-empty-file "
               "--add-index=tracking_id "
               "--add-index=control_name "
               "--add-index=test_id")

    # Jethro - load tables of sample specific cuffdiff fpkm values into csvdb
    # IMS: First read in lookup table for CuffDiff/Pipeline sample name
    # conversion
    inf = iotools.open_file(os.path.join(indir, "read_groups.info.gz"))
    inf.readline()
    sample_lookup = {}

    for line in inf:
        line = line.split("\t")
        our_sample_name = iotools.snip(line[0])
        our_sample_name = re.sub("-", "_", our_sample_name)
        cuffdiff_sample_name = "%s_%s" % (line[1], line[2])
        sample_lookup[cuffdiff_sample_name] = our_sample_name

    inf.close()

    for fn, level in (("cds.read_group_tracking.gz", "cds"),
                      ("genes.read_group_tracking.gz", "gene"),
                      ("isoforms.read_group_tracking.gz", "isoform"),
                      ("tss_groups.read_group_tracking.gz", "tss")):

        tablename = prefix + "_" + level + "sample_fpkms"

        tmpf = P.get_temp_filename(".")
        inf = iotools.open_file(os.path.join(indir, fn)).readlines()
        outf = iotools.open_file(tmpf, "w")

        samples = []
        genes = {}

        is_first = True
        for line in inf:

            if is_first:
                is_first = False
                continue

            line = line.split()
            gene_id = line[0]
            condition = line[1]
            replicate = line[2]
            fpkm = line[6]
            status = line[8]

            sample_id = condition + "_" + replicate

            if sample_id not in samples:
                samples.append(sample_id)

            # IMS: The following block keeps getting its indenting messed
            # up. It is not part of the 'if sample_id not in samples' block
            # please make sure it does not get made part of it
            if gene_id not in genes:
                genes[gene_id] = {}
                genes[gene_id][sample_id] = fpkm
            else:
                if sample_id in genes[gene_id]:
                    raise ValueError(
                        'sample_id %s appears twice in file for gene_id %s'
                        % (sample_id, gene_id))
                else:
                    if status != "OK":
                        genes[gene_id][sample_id] = status
                    else:
                        genes[gene_id][sample_id] = fpkm

        samples = sorted(samples)

        # IMS - CDS files might be empty if not cds has been
        # calculated for the genes in the long term need to add CDS
        # annotation to denovo predicted genesets in meantime just
        # skip if cds tracking file is empty

        if len(samples) == 0:
            continue

        headers = "gene_id\t" + "\t".join([sample_lookup[x] for x in samples])
        outf.write(headers + "\n")

        for gene in genes.iterkeys():
            outf.write(gene + "\t")
            s = 0
            while x < len(samples) - 1:
                outf.write(genes[gene][samples[s]] + "\t")
                s += 1

            # IMS: Please be careful with this line. It keeps getting moved
            # into the above while block where it does not belong
            outf.write(genes[gene][samples[len(samples) - 1]] + "\n")

        outf.close()

        P.load(tmpf,
               outfile,
               tablename=tablename,
               options="--allow-empty-file "
               " --add-index=gene_id")

        os.unlink(tmpf)

    # build convenience table with tracks
    tablename = prefix + "_isoform_levels"
    tracks = Database.getColumnNames(dbhandle, tablename)
    tracks = [x[:-len("_FPKM")] for x in tracks if x.endswith("_FPKM")]

    tmpfile = P.get_temp_file(dir=".")
    tmpfile.write("track\n")
    tmpfile.write("\n".join(tracks) + "\n")
    tmpfile.close()

    P.load(tmpfile.name, outfile)
    os.unlink(tmpfile.name)


def parseCuffdiff(infile, min_fpkm=1.0):
    '''parse a cuffdiff .diff output file.

    This method takes cuffdiff output and converts the results into a
    standardized table.

    Arguments
    ---------
    infile : string
        Input filename, output from cuffdiff
    min_fpkm : float
        Minimum fpkm. Genes with an fpkm lower than this will
        be set to status `NOCALL`.
    '''

    CuffdiffResult = collections.namedtuple(
        "CuffdiffResult",
        "test_id gene_id gene  locus   sample_1 sample_2  "
        "status  value_1 value_2 l2fold  "
        "test_stat p_value q_value significant ")

    results = []

    for line in iotools.open_file(infile):
        if line.startswith("test_id"):
            continue
        data = CuffdiffResult._make(line[:-1].split("\t"))
        status = data.status
        significant = [0, 1][data.significant == "yes"]
        if status == "OK" and (float(data.value_1) < min_fpkm or
                               float(data.value_2) < min_fpkm):
            status = "NOCALL"

        try:
            fold = math.pow(2.0, float(data.l2fold))
        except OverflowError:
            fold = "na"

        results.append(Expression.GeneExpressionResult._make((
            data.test_id,
            data.sample_1,
            data.value_1,
            0,
            data.sample_2,
            data.value_2,
            0,
            data.p_value,
            data.q_value,
            data.l2fold,
            fold,
            data.l2fold,
            significant,
            status)))

    return results


def runCuffdiff(bamfiles,
                design_file,
                geneset_file,
                outfile,
                cuffdiff_options="",
                job_threads=4,
                job_memory="4G",
                fdr=0.1,
                mask_file=None):
    '''estimate differential expression using cuffdiff.

    Replicates within each track are grouped.

    Arguments
    ---------
    bamfiles : list
        List of filenames in :term:`bam` format.
    designfile : string
        Filename with experimental design in :term:`tsv` format.
    geneset_file : string
        Filename with geneset of interest in :term:`gtf format.
    outfile : string
        Output filename. The output is :term:`tsv` formatted.
    cuffdiff_options : string
        Options to pass on to cuffdiff
    job_threads : int
        Number of threads to use.
    job_memory : string
        Memory to reserve.
    fdr : float
        FDR threshold to apply.
    mask_file : string
        If given, ignore genes overlapping gene models in
        this :term:`gtf` formatted file.
    '''

    design = Expression.readDesignFile(design_file)

    outdir = outfile + ".dir"
    try:
        os.mkdir(outdir)
    except OSError:
        pass

    # replicates are separated by ","
    reps = collections.defaultdict(list)
    for bamfile in bamfiles:
        groups = collections.defaultdict()
        # .accepted.bam kept for legacy reasons (see rnaseq pipeline)
        track = P.snip(os.path.basename(bamfile), ".bam", ".accepted.bam")
        if track not in design:
            E.warn("bamfile '%s' not part of design - skipped" % bamfile)
            continue

        d = design[track]
        if not d.include:
            continue
        reps[d.group].append(bamfile)

    groups = sorted(reps.keys())
    labels = ",".join(groups)
    reps = "   ".join([",".join(reps[group]) for group in groups])

    # Nick - add mask gtf to not assess rRNA and ChrM
    extra_options = []

    if mask_file:
        extra_options.append(" -M %s" % os.path.abspath(mask_file))

    extra_options = " ".join(extra_options)

    # AH: removed log messages about BAM record error
    # These cause logfiles to grow several Gigs and are
    # frequent for BAM files not created by tophat.
    # Error is:
    # BAM record error: found spliced alignment without XS attribute
    # AH: compress output in outdir
    job_memory = "7G"
    statement = '''date > %(outfile)s.log;
    hostname >> %(outfile)s.log;
    cuffdiff --output-dir %(outdir)s
             --verbose
             --num-threads %(job_threads)i
             --labels %(labels)s
             --FDR %(fdr)f
             %(extra_options)s
             %(cuffdiff_options)s
             <(gunzip < %(geneset_file)s )
             %(reps)s
    2>&1
    | grep -v 'BAM record error'
    >> %(outfile)s.log;
    gzip -f %(outdir)s/*;
    date >> %(outfile)s.log;
    '''
    P.run(statement)

    results = parseCuffdiff(os.path.join(outdir, "gene_exp.diff.gz"))

    Expression.writeExpressionResults(outfile, results)


# UTR estimation
Utr = collections.namedtuple("Utr", "old new max status")


def buildUTRExtension(infile, outfile):
    '''build new utrs by building and fitting an HMM
    to reads upstream and downstream of known genes.

    Known problems

    * the size of the extension is limited by the window size

    * introns within UTRs are ignored.

    * UTR extension might be underestimated for highly expressed genes
      as relative read counts drop off quickly, even though there is
      a good amount of reads still present in the UTR.

    The model

    The model is a three-state model::

        UTR --|--> notUTR --|--> otherTranscript --|
          ^---|      ^------|              ^-------|
                     ^-----------------------------|

    The chain starts in UTR and ends in notUTr or otherTranscript.

    The otherTranscript state models peaks of within the upstream/
    downstream region of a gene. These peaks might correspond to
    additional exons or unknown transcripts. Without this state,
    the UTR might be artificially extend to include these peaks.

    Emissions are modelled with beta distributions. These
    distributions permit both bimodal (UTR) and unimodal (notUTR)
    distribution of counts.

    Parameter estimation

    Parameters are derived from known UTRs within full length
    territories.

    Transitions and emissions for the otherTranscript state
    are set heuristically:

       * low probabibily for remaining in state "otherTranscript".
           * these transcripts should be short.

       * emissions biased towards high counts - only strong signals
           will be considered.

       * these could be estimated from known UTRs, but I am worried
           UTR extensions then will be diluted.


    Alternatives

    The method could be improved.

    * base level resolution?
       * longer chains result in more data and longer running times.
       * the averaging in windows smoothes the data, which might have
         a beneficial effect.

    * raw counts instead of scaled counts?
       * better model, as highly expressed genes should give more
         confident predictions.

    Arguments
    ---------
    infile : string
        Output of :func:`buildGeneLevelReadExtension`
    outfile : string
        Output filename

    '''

    # the bin size , see gtf2table - can be cleaned from column names
    # or better set as options in .yml file
    binsize = 100
    territory_size = 15000

    # read gene coordinates
    geneinfos = {}
    for x in CSV.DictReader(iotools.open_file(infile), dialect='excel-tab'):
        contig, strand, start, end = x['contig'], x[
            'strand'], int(x['start']), int(x['end'])
        geneinfos[x['gene_id']] = (contig, strand,
                                   start, end)

    infiles = [infile + ".readextension_upstream_sense.tsv.gz",
               infile + ".readextension_downstream_sense.tsv.gz"]

    outdir = os.path.join(PARAMS["exportdir"], "utr_extension")

    R('''suppressMessages(library(RColorBrewer))''')
    R('''suppressMessages(library(MASS))''')
    R('''suppressMessages(library(HiddenMarkov))''')

    # for upstream, downstream
    upstream_utrs, downstream_utrs = {}, {}

    all_genes = set()

    for filename, new_utrs in zip(infiles, (upstream_utrs, downstream_utrs)):

        E.info("processing %s" % filename)

        parts = os.path.basename(filename).split(".")

        data = R(
            '''data = read.table(gzfile( "%(filename)s"), header=TRUE,
            fill=TRUE, row.names=1)''' % locals())

        ##########################################
        ##########################################
        ##########################################
        # estimation
        ##########################################
        # take only those with a 'complete' territory
        R('''d = data[-which( apply( data,1,function(x)any(is.na(x)))),]''')
        # save UTR
        R('''utrs = d$utr''')
        # remove length and utr column
        R('''d = d[-c(1,2)]''')
        # remove those which are completely empty, logtransform or scale data
        # and export
        R('''lraw = log10(
        d[-which(apply(d, 1, function(x)all(x==0))),] + 1)''')

        utrs = R('''utrs = utrs[-which( apply(d,1,function(x)all(x==0)))]''')
        scaled = R(
            '''lscaled = t(scale(t(lraw), center=FALSE,
            scale=apply(lraw,1,max)))''')
        exons = R('''lraw[,1]''')

        #######################################################
        #######################################################
        #######################################################
        # do the estimation:
        E.debug("estimation: utrs=%i, exons=%i, vals=%i, dim=%s" %
                (len(utrs), len(exons), len(scaled), R.dim(scaled)))
        # counts within and outside UTRs
        within_utr, outside_utr, otherTranscript = [], [], []
        # number of transitions between utrs
        transitions = np.zeros((3, 3), np.int)

        for x in xrange(len(utrs)):
            utr, exon = utrs[x], exons[x]

            # only consider genes with expression coverage
            # note: expression level is logscaled here, 10^1 = 10
            if exon < 0.1:
                continue

            # first row is column names, so x + 1
            values = list(scaled.rx(x + 1, True))

            utr_bins = utr // binsize
            nonutr_bins = (territory_size - utr) // binsize

            # build transition matrix
            transitions[0][0] += utr_bins
            transitions[0][1] += 1
            transitions[1][1] += nonutr_bins

            outside_utr.extend([x for x in values[utr_bins:] if x <= 0.5])

            # ignore exon and zero counts
            within_utr.extend([x for x in values[1:utr_bins] if x > 0.1])

            # add only high counts to otherTranscript emissions
            otherTranscript.extend([x for x in values[utr_bins:] if x > 0.5])

        # estimation for
        # 5% chance of transiting to otherTranscript
        transitions[1][2] = transitions[1][1] * 0.05
        # 10% chance of remaining in otherTranscript
        transitions[2][1] = 900
        transitions[2][2] = 100

        E.info("counting: (n,mean): within utr=%i,%f, "
               "outside utr=%i,%f, otherTranscript=%i,%f" %
               (len(within_utr), np.mean(within_utr),
                len(outside_utr), np.mean(outside_utr),
                len(otherTranscript), np.mean(otherTranscript)))

        ro.globalenv['transitions'] = R.matrix(transitions, nrow=3, ncol=3)
        R('''transitions = transitions / rowSums( transitions )''')
        ro.globalenv['within_utr'] = ro.FloatVector(within_utr[:10000])
        ro.globalenv['outside_utr'] = ro.FloatVector(outside_utr[:10000])
        ro.globalenv['otherTranscript'] = ro.FloatVector(
            otherTranscript[:10000])

        # estimate beta distribution parameters
        R('''doFit = function( data ) {
                   data[data == 0] = data[data == 0] + 0.001
                   data[data == 1] = data[data == 1] - 0.001
                   f = fitdistr( data, dbeta, list( shape1=0.5, shape2=0.5 ) )
                   return (f) }''')

        fit_within_utr = R(
            '''fit_within_utr = suppressMessages(doFit( within_utr))''')
        fit_outside_utr = R(
            '''fit_outside_utr = suppressMessages(doFit( outside_utr))''')
        fit_other = R(
            '''fit_otherTranscript = suppressMessages(
            doFit(otherTranscript))''')

        within_a, within_b = list(fit_within_utr.rx("estimate"))[0]
        outside_a, outside_b = list(fit_outside_utr.rx("estimate"))[0]
        other_a, other_b = list(fit_other.rx("estimate"))[0]

        E.info("beta estimates: within_utr=%f,%f outside=%f,%f, other=%f,%f" %
               (within_a, within_b, outside_a, outside_b, other_a, other_b))

        fn = ".".join((parts[0], parts[4], "fit", "png"))
        outfilename = os.path.join(outdir, fn)
        R.png(outfilename, height=1000, width=1000)

        R('''par(mfrow=c(3,1))''')
        R('''x=seq(0,1,0.02)''')
        R('''hist( within_utr, 50, col=rgb( 0,0,1,0.2) )''')
        R('''par(new=TRUE)''')
        R('''plot(x, dbeta(x, fit_within_utr$estimate['shape1'],
        fit_within_utr$estimate['shape2']), type='l', col='blue')''')

        R('''hist( outside_utr, 50, col=rgb( 1,0,0,0.2 ) )''')
        R('''par(new=TRUE)''')
        R('''plot( x, dbeta( x, fit_outside_utr$estimate['shape1'],
        fit_outside_utr$estimate['shape2']), type='l', col='red')''')

        R('''hist( otherTranscript, 50, col=rgb( 0,1,0,0.2 ) )''')
        R('''par(new=TRUE)''')
        R('''plot( x, dbeta( x, fit_otherTranscript$estimate['shape1'],
        fit_otherTranscript$estimate['shape2']), type='l', col='green')''')
        R['dev.off']()

        #####################################################
        #####################################################
        #####################################################
        # build hmm
        # state 1 = UTR
        # state 2 = notUTR
        # state 3 = other transcript
        p = R('''betaparams = list( shape1=c(fit_within_utr$estimate['shape1'],
        fit_outside_utr$estimate['shape1'],
        fit_otherTranscript$estimate['shape1']),
        shape2=c(fit_within_utr$estimate['shape2'],
        fit_outside_utr$estimate['shape2'],
        fit_otherTranscript$estimate['shape2'])) ''')
        R('''hmm = dthmm(NULL, transitions, c(1,0,0), "beta", betaparams )''')

        E.info("fitting starts")
        #####################################################
        #####################################################
        #####################################################
        # fit to every sequence
        genes = R('''rownames(data)''')
        all_genes.update(set(genes))
        utrs = R('''data$utr''')
        exons = R('''data$exon''')
        nseqs = len(utrs)

        counter = E.Counter()

        for idx in xrange(len(utrs)):

            gene_id = genes[idx]

            old_utr = utrs[idx]

            if idx % 100 == 0:
                E.debug("processing gene %i/%i" % (idx, len(utrs)))

            counter.input += 1

            # do not predict if terminal exon not expressed
            if exons[idx] < 1:
                counter.skipped_notexpressed += 1
                new_utrs[gene_id] = Utr._make(
                    (old_utr, None, None, "notexpressed"))
                continue

            R('''obs = data[%i,][-c(1,2)]''' % (idx + 1))
            # remove na
            obs = R('''obs = obs[!is.na(obs)]''')
            if len(obs) <= 1 or max(obs) == 0:
                new_utrs[gene_id] = Utr._make(
                    (old_utr, None, None, "no observations"))
                continue

            # normalize
            R('''obs = obs / max(obs)''')
            # add small epsilon to 0 and 1 values
            R('''obs[obs==0] = obs[obs==0] + 0.001 ''')
            R('''obs[obs==1] = obs[obs==1] - 0.001 ''')
            R('''hmm$x = obs''')

            states = None
            try:
                states = list(R('''states = Viterbi( hmm )'''))
            except ri.RRuntimeError as msg:
                counter.skipped_error += 1
                new_utrs[gene_id] = Utr._make((old_utr, None, None, "fail"))
                continue

            max_utr = binsize * (len(states) - 1)

            # subtract 1 for last exon
            try:
                new_utr = binsize * (states.index(2) - 1)
                new_utrs[gene_id] = Utr._make(
                    (old_utr, new_utr, max_utr, "ok"))
                counter.success += 1
            except ValueError:
                new_utrs[gene_id] = Utr._make(
                    (old_utr, max_utr, max_utr, "max"))
                counter.maxutr += 1

    E.info("fitting: %s" % str(counter))

    outf = iotools.open_file(outfile, "w")

    outf.write("\t".join(
        ["gene_id", "contig", "strand", "status5", "status3"] +
        ["%s_%s_%s" % (x, y, z) for x, y, z in itertools.product(
            ("old", "new", "max"),
            ("5utr", "3utr"),
            ("length", "start", "end"))]) + "\n")

    def _write(coords, strand):

        start5, end5, start3, end3 = coords
        if strand == "-":
            start5, end5, start3, end3 = start3, end3, start5, end5

        if start5 is None:
            start5, end5, l5 = "", "", ""
        else:
            l5 = end5 - start5

        if start3 is None:
            start3, end3, l3 = "", "", ""
        else:
            l3 = end3 - start3

        return "\t".join(map(str, (l5, start5, end5,
                                   l3, start3, end3)))

    def _buildCoords(upstream, downstream, start, end):

        if upstream:
            start5, end5 = start - upstream, start
        else:
            start5, end5 = None, None
        if downstream:
            start3, end3 = end, end + downstream
        else:
            start3, end3 = None, None

        return start5, end5, start3, end3

    for gene_id in all_genes:

        contig, strand, start, end = geneinfos[gene_id]

        outf.write("%s\t%s\t%s" % (gene_id, contig, strand))

        if gene_id in upstream_utrs:
            upstream = upstream_utrs[gene_id]
        else:
            upstream = Utr._make((None, None, None, "missing"))
        if gene_id in downstream_utrs:
            downstream = downstream_utrs[gene_id]
        else:
            downstream = Utr._make((None, None, None, "missing"))

        if strand == "-":
            upstream, downstream = downstream, upstream

        # output prediction status
        outf.write("\t%s\t%s" % (upstream.status, downstream.status))

        # build upstream/downstream coordinates
        old_coordinates = _buildCoords(
            upstream.old, downstream.old, start, end)
        new_coordinates = _buildCoords(
            upstream.new, downstream.new, start, end)

        # reconciled = take maximum extension of UTR
        max_coordinates = []
        # note that None counts as 0 in min/max.
        for i, d in enumerate(zip(old_coordinates, new_coordinates)):
            if i % 2 == 0:
                v = [z for z in d if z is not None]
                if v:
                    max_coordinates.append(min(v))
                else:
                    max_coordinates.append(None)
            else:
                max_coordinates.append(max(d))

        # convert to 5'/3' coordinates
        outf.write("\t%s\t%s\t%s\n" % (_write(old_coordinates, strand),
                                       _write(new_coordinates, strand),
                                       _write(max_coordinates, strand)))

    outf.close()


def plotGeneLevelReadExtension(infile, outfile):
    '''plot read density of reads extending beyond last exon.

    Arguments
    ---------
    infile : string
    outfile : string

    '''

    infiles = glob.glob(infile + ".*.tsv.gz")

    R('''suppressMessages(library(RColorBrewer))''')
    R('''suppressMessages(library(MASS))''')
    R('''suppressMessages(library(HiddenMarkov))''')

    # the bin size , see gtf2table - could be cleaned from column names
    binsize = 100
    territory_size = 15000

    for filename in infiles:

        E.info("processing %s" % filename)

        parts = os.path.basename(filename).split(".")

        data = R(
            '''data = read.table(gzfile("%(filename)s"),
            header=TRUE, fill=TRUE, row.names=1)''' % locals())

        # estimation
        # take only those with a 'complete' territory
        R('''d = data[-which( apply( data,1,function(x)any(is.na(x)))),]''')
        # save UTR
        R('''utrs = d$utr''')
        # remove length and utr column
        R('''d = d[-c(1,2)]''')
        # remove those which are completely empty, logtransform or scale data
        # and export
        R('''lraw = log10(d[-which( apply(d,1,function(x)all(x==0))),] + 1)''')

        utrs = R('''utrs = utrs[-which( apply(d,1,function(x)all(x==0)))]''')
        scaled = R(
            '''lscaled = t(scale(t(lraw), center=FALSE,
            scale=apply(lraw,1,max)))''')
        exons = R('''lraw[,1]''')

        if len(utrs) == 0:
            E.warn("no data for %s" % filename)
            continue

        #######################################################
        #######################################################
        #######################################################
        R('''myplot = function( reads, utrs, ... ) {
           oreads = t(data.matrix( reads )[order(utrs), ] )
           outrs = utrs[order(utrs)]
           image( 1:nrow(oreads), 1:ncol(oreads), oreads ,
                  xlab = "", ylab = "",
                  col=brewer.pal(9,"Greens"),
                  axes=FALSE)
           # axis(BELOW<-1, at=1:nrow(oreads), labels=rownames(oreads),
           # cex.axis=0.7)
           par(new=TRUE)
           plot( outrs, 1:length(outrs), yaxs="i", xaxs="i",
                 ylab="genes", xlab="len(utr) / bp",
                 type="S",
                 xlim=c(0,nrow(oreads)*%(binsize)i))
        }''' % locals())

        fn = ".".join((outfile, parts[0], parts[4], "raw", "png"))

        R.png(fn, height=2000, width=1000)
        R('''myplot(lraw, utrs)''')
        R['dev.off']()

        # plot scaled data
        fn = ".".join((outfile, parts[0], parts[4], "scaled", "png"))

        R.png(fn, height=2000, width=1000)
        R('''myplot(lscaled, utrs)''')
        R['dev.off']()

    iotools.touch_file(outfile)
