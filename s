[33mcommit 4d66f049115e68cb8722504d754aea5d1b5fabeb[m[33m ([m[1;36mHEAD -> [m[1;32m73-options-sample-composition-graphs[m[33m)[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Jun 21 12:20:24 2021 +0200

    methods reorganization inside PlotTaxonomyCounts

[33mcommit d48d5254aae8cc3efa27853f66428621b24f8725[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu May 6 14:57:55 2021 +0200

    update documentation link

[33mcommit 08ca8e477f663d2e641f224d00de956311953589[m
Merge: fddf1f3 db7acc8
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu May 6 14:54:06 2021 +0200

    Merge pull request #72 from motleystate/71-quick-taxonomy-viz
    
    Built in visualization for taxonomy counts

[33mcommit db7acc872bde4365450d394dd40314a8194a0aaa[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu May 6 13:58:57 2021 +0200

    fix automatic Y axis label
    
    update docstring

[33mcommit 7b2787d3bb789aa8d72d06ad00811ba333522b58[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed May 5 14:14:43 2021 +0200

    optinal clustering and possibility to force order

[33mcommit bbb818e965dcd41f8c99e25bec56f1ed8f085b70[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed May 5 11:32:58 2021 +0200

    refactor and start testing private methods for PlotTaxonomyCounts

[33mcommit 76828e55020f69789c0dbd01c20e7e7c14ba542f[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue May 4 14:29:14 2021 +0200

    merge plotting options

[33mcommit 830fdb1a09cb5b4a4c29af2d495620cfed764788[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue May 4 13:38:31 2021 +0200

    code gardening

[33mcommit 83005f9be773d3e20a5d9e9e0f972666f31fb21e[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue May 4 13:37:46 2021 +0200

    display mean of taxa among samples

[33mcommit 02bc68905aac2f3c66f26a6ee9b8a1d02fd63630[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue May 4 11:33:23 2021 +0200

    fix trace order for legend

[33mcommit b0cff15c3c3683a5c1aa88e52c21da4809df2f7f[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue May 4 10:22:45 2021 +0200

    fix plotting options issue for plot_sample_composition_most_abundant_taxa
    
    plotting_options for plot_most_abundant_taxa

[33mcommit 3bee5cb72bc50e845ed3b25d29dc682a8e41f5da[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon May 3 16:51:11 2021 +0200

    fix flake8

[33mcommit 8fd5be7e42d1d5527cb320c3086dbe9e35e664ab[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon May 3 16:49:06 2021 +0200

    automatic colors on name for bar charts

[33mcommit 394c6519c2878043cae33e19990f2910326aa333[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Fri Apr 30 17:46:24 2021 +0200

    start method to plot sample compositions

[33mcommit 87fd8aefe652bc76b5b137c6422ab0c4dc6e552a[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Apr 29 15:33:53 2021 +0200

    option to MatrixBarGraph for colors of group (index)

[33mcommit 351f0ca33d7ab2c0aad172010c7e41f2b5673ac7[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Apr 29 15:24:36 2021 +0200

    add BarGraph to display matrix data

[33mcommit 2a7097f120b621810aad959d17972add48071191[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Apr 29 11:23:57 2021 +0200

    rename and make loading data private method

[33mcommit 1e2c6095621306b9f319e29cc3b40617fd5a6727[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Apr 29 10:28:05 2021 +0200

    create PlotTaxonomyCounts and general way to plot from parser

[33mcommit 653c7cda60a4a14248a3efe1daf6c0ed81088cef[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed Apr 28 12:44:51 2021 +0200

    visualize top taxa for taxonomy counts
    
    run black on new code

[33mcommit 9a4dbb1670fc42576031e7c6562d9704589f3954[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue Apr 27 17:46:22 2021 +0200

    add TaxonomyMeanFiltering to filtering imports

[33mcommit 5bb81073773e0976f4848bd62b05201fc4da4a5f[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue Apr 27 17:38:44 2021 +0200

    add filtering on mean among samples

[33mcommit e2720a1aacdebfafb3e64390d0444a7be102182b[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Tue Apr 27 11:43:28 2021 +0200

    reorganize taxonomy parsers and create common base

[33mcommit fddf1f3d810291ea0438b7810af4cd525ced4698[m
Merge: 709c690 1647d9c
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Wed Apr 7 12:33:45 2021 +0200

    Merge pull request #66 from motleystate/adapting-metaphlan3-parser-for-relative-abundance
    
    adaptation of metaphlan3 parser

[33mcommit 1647d9c32ae04754a59cfc38f2bd6b6011763661[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed Apr 7 12:27:18 2021 +0200

    remove scikit-learn version

[33mcommit a160bed24ccbf67b0671973317c98d7a350830a2[m
Merge: 823769d 709c690
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed Apr 7 12:27:04 2021 +0200

    Merge branch 'master' into adapting-metaphlan3-parser-for-relative-abundance

[33mcommit 709c690c3b9d442cde25c65c815ade6b6bd2f556[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Apr 1 12:14:17 2021 +0200

    Chao1Index

[33mcommit 057622f013696e955c93e9f4b40b605ca935436b[m
Merge: 4c4d6c1 941cde8
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Tue Feb 16 10:29:37 2021 +0100

    Merge pull request #70 from motleystate/small-fixes-and-additional-options
    
    Small fixes and additional options

[33mcommit 941cde85c2e835b3bf74b2260862d416126b90a5[m
Merge: 6754b29 a58e9c7
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Feb 15 18:42:32 2021 +0100

    Merge branch 'small-fixes-and-additional-options' of https://github.com/motleystate/moonstone into small-fixes-and-additional-options

[33mcommit 6754b29a0205f199e334b4265028ffc392c5a591[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Feb 15 18:42:09 2021 +0100

    test np.nan fix

[33mcommit 823769d5540ac9a3e8b542a418785f01b49d2698[m
Merge: 23c390f 4c4d6c1
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon Feb 15 11:20:14 2021 +0100

    Merge branch 'master' into adapting-metaphlan3-parser-for-relative-abundance

[33mcommit a58e9c7d6e40255f1c95202114d7ce1f7e3a9cd1[m
Merge: b97d5b8 4c4d6c1
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon Feb 15 10:43:17 2021 +0100

    Merge branch 'master' into small-fixes-and-additional-options

[33mcommit 4c4d6c18bf7d1dc888f7c0c35168d91ec90e2979[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon Feb 15 10:13:43 2021 +0100

    fix: fix CI by fixing hdmedians (0.14.1 leads to issues)

[33mcommit df7b9f7d3f63c7a3fbd86bcaba443b3ea6740dee[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon Feb 15 09:55:29 2021 +0100

    fix: install right version of numpy in CI

[33mcommit 7dc7a145d06f4458838657a2de8bc24b4fd4ae30[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon Feb 15 09:45:45 2021 +0100

    fix: fix versions for tests in CI

[33mcommit 3433abfe06213f2ddc8e31c40702e2ca3970fdc6[m
Merge: 83163f2 4a35cfb
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Fri Feb 12 21:08:37 2021 +0100

    Merge pull request #69 from motleystate/add-counts-bar-plots
    
    add group counts for group Box

[33mcommit b97d5b87977d4aa1b0b3babf1b1df4ca0331ad10[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Feb 9 12:00:09 2021 +0100

    fix correction pval when NaN + option png output

[33mcommit 83163f26d9205135bc060e505201579795b31552[m[33m ([m[1;33mtag: v0.5.0[m[33m)[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Jan 28 14:41:59 2021 +0100

    bump to version 0.5.0

[33mcommit 4a35cfb933a854c601821d857764956261b4abd0[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Mon Jan 25 16:52:37 2021 +0100

    add group counts for group Box

[33mcommit b0610bbb6bfe6bab2e14c91408f5606fdcc5c4fe[m
Merge: 3e88649 61a131f
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Thu Jan 21 12:05:30 2021 +0100

    Merge pull request #68 from motleystate/67-_get_grouped_df-np.nan-bug
    
    fixes bug with np.nan in BetaDiversity

[33mcommit 61a131f1b17f6a07ffa12dfdb306ab630c86490b[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Jan 21 11:58:45 2021 +0100

    fixes bug

[33mcommit 23c390ff121e8dc65ea1e6d75f4666d909f21742[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Jan 19 15:13:34 2021 +0100

    last corrections (hopefully)

[33mcommit 3f4442a698e373fe2858eb8263f639a089170fb7[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Jan 14 19:05:22 2021 +0100

    in metaphlan2 + rank_level dilemma

[33mcommit b9252d60fae503ab7d5f89e80e903748f3ddd7b8[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Jan 13 19:02:34 2021 +0100

    corrections + creation of BaseMetaphlanParser

[33mcommit 10b3c8f338c2c98238928d1d831a4c26f347ce64[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Jan 11 14:40:55 2021 +0100

    adapting metaphlan3

[33mcommit 3e886497b5726211ed8b0cc74ba9e7ceec9becb5[m
Merge: aa19607 11a4007
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Thu Jan 7 11:23:15 2021 +0100

    Merge pull request #62 from motleystate/visualize-metadata-cat
    
    Visualize metadata cat

[33mcommit 11a40072981b7b4e0020c5e1827291128446f3a1[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Jan 7 11:11:56 2021 +0100

    add output file and title to graph

[33mcommit bc306582e806352672c271b16c3336dd594f340d[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed Jan 6 16:16:49 2021 +0100

    code gardening

[33mcommit 8063c46652bf2068afd82c9e603bee11c0266813[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed Jan 6 12:40:46 2021 +0100

    add tests for metadata parser

[33mcommit fd1651ec70dcbb3c9f9a713e1e76e806b37f3c95[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Wed Jan 6 11:52:25 2021 +0100

    update docs

[33mcommit 320198aa14690268eea8fc5c7df853c457e19d71[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Dec 31 16:48:15 2020 +0100

    remove font style

[33mcommit 2d1d99f65740934f202b02a7a3bff3e131bb314d[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@posteo.net>
Date:   Thu Dec 31 16:14:03 2020 +0100

    add method to visualize categories of metadata

[33mcommit aa196078c9cec27b54ee117939697d698be4829f[m
Merge: f84954b b1266c7
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Thu Dec 10 15:28:53 2020 +0100

    Merge pull request #61 from motleystate/59-solving-error-NamesFiltering
    
    Correction of error in NamesFiltering due to pandas deprecation

[33mcommit b1266c75d813bb21e735d31db7e9b9ed8d97ebd2[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Dec 10 15:18:05 2020 +0100

    only 1 warning msg

[33mcommit c3e97b748a1decc9f07a23557af9a6af64c31f89[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Dec 10 12:16:09 2020 +0100

    error resolution (hopefully) + tests

[33mcommit f84954bf512ad91c7c261c656d06e8ab6899bd05[m
Merge: 5f29984 2b106ae
Author: skennedy <seanpkennedy@gmail.com>
Date:   Thu Dec 10 10:35:21 2020 +0100

    Merge pull request #60 from motleystate/metaphlan3-parser
    
    start parser for metaphlan3

[33mcommit 2b106ae83e6cd006a88075525bcf882dde4db5a7[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Wed Dec 9 18:39:05 2020 +0100

    update documentation

[33mcommit ed05c3076112d0a0bf6ca40653d279e558a54e8c[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Wed Dec 9 18:36:54 2020 +0100

    start parser for metaphlan3

[33mcommit 5f2998499e246e3ee661115cd04bbe3727d3def3[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Nov 27 15:36:39 2020 +0100

    Some additional optimization: total reads and an x-axis title.
    x-axis in now log-scale

[33mcommit ccf145e9e0e2c3cb17b6bcf968ea98b48a6d0140[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Nov 27 10:26:32 2020 +0100

    Finished with `make_annotations` function.

[33mcommit 083c884cfc14354a7b03223c5f6429ae36871c8c[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Nov 19 11:30:09 2020 +0100

    WIP: Started `plat_placeholder` module for later integration with
    other plotting functions.
    
    Collects read info and generates annotation, formats boxplot.

[33mcommit 390d0f81f672e2b6d1c8cbd81f23b57d66237995[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Nov 13 14:41:20 2020 +0100

    Added `plot_reads` function.
    Still need to add the actual plotting part.

[33mcommit 4cf1144b66d0b0097e003cf3c4d322bd3643b34e[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Nov 13 11:27:34 2020 +0100

    Added further explanatory details to `read_info` function.

[33mcommit d2fd6e121b22044e724e43a47f79b6adce4215db[m
Merge: 909b2f0 6fe4926
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Thu Nov 12 14:32:35 2020 +0100

    Merge pull request #57 from motleystate/pvalue-correction
    
    pvalue correction in analyse groups

[33mcommit 6fe49266d1ad354f969aea9af4b29b6801d2eb5e[m
Merge: 06d7086 909b2f0
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Nov 10 12:41:51 2020 +0100

    Merge branch 'master' into pvalue-correction

[33mcommit 06d70863307d11a59b40e5033b16743ae802732a[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Nov 10 11:51:02 2020 +0100

    bug corrections

[33mcommit 909b2f0558257b1ccd4d9880a81bd022d7ac2042[m
Merge: 74db716 60e257e
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Tue Nov 10 11:42:26 2020 +0100

    Merge pull request #56 from motleystate/55-heatmap
    
    add HeatmapGraph

[33mcommit 85fcb475b71c7caafceab8a2cc60d7b9280e18c3[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Nov 10 11:33:03 2020 +0100

    tests + slight refactoring

[33mcommit 74db716d3084f1be01782a66bd5bb90e247af677[m
Merge: 7db21f5 34ed1dd
Author: skennedy <seanpkennedy@gmail.com>
Date:   Tue Nov 10 11:23:34 2020 +0100

    Merge pull request #54 from motleystate/52_MPsupport
    
    52 m psupport

[33mcommit dab874fb118c688eb2420299fb13d48f2e1ff1c5[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Nov 9 12:04:00 2020 +0100

    pvalue correction in analyse groups

[33mcommit 60e257ec7048df6874609c1aed88b7c60b2c3b7a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Nov 9 11:19:43 2020 +0100

    handle output file for pvalue heatmaps

[33mcommit 24cffd8a72e1c574f9e21d3381c281082bfc4691[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Nov 9 11:16:47 2020 +0100

    use heatmap to visualize pval from diversity

[33mcommit 31900f614f05434ef0e265177365abbedf18dc90[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Nov 9 11:00:29 2020 +0100

    add HeatmapGraph

[33mcommit 34ed1dd2b7654b1dd911b5de569b90f08ac08051[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Nov 5 14:42:07 2020 +0100

    Reset to `origin/master`

[33mcommit 5317b83db1c2b1238f2c3e902fefa0894b7bacc8[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Nov 5 14:20:06 2020 +0100

    Changed handling of uncompressed files which normally yielded a NONE type
    from 'filetype`. This is now explicitly set to 'Uncompressed/FASTQ' to match
    `downsize_dir` functions.
    
    Corrected line in `count_starting_reads` that redundantly called `filetype`
    
    Incorporated explicit filetype in the `downsize_pair` function.
    Corrected calling the wrong function here.

[33mcommit 23db139fb8f87e1340dc91788b0c119112dbf10e[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Nov 5 14:15:15 2020 +0100

    Added module to handle setting parameters for each worker in the
    multiprocessing pool.
    
    Clarified some variables naming.

[33mcommit d269314f8f0d2a38b68f34acdfc6a3a189aea28d[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Nov 5 14:08:31 2020 +0100

    rename to match `read_downsize.py` function.

[33mcommit 7db21f5ebf411d7453546cfce2c7109adbe50b42[m
Merge: e81f93d 5c3364c
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Thu Nov 5 11:01:22 2020 +0100

    Merge pull request #49 from motleystate/48-beta-diversity
    
    Beta Diversity

[33mcommit 5c3364c7c1613bd84199869f34532d3d63c17d12[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 5 10:40:15 2020 +0100

    add meta to remind the stats test used

[33mcommit 3218a4e7d3ebfe7674b083f33ecf28e602400fcb[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 5 10:26:17 2020 +0100

    remove duplicated code

[33mcommit 292b1b4a67ba69f7b7641641029c9e9b8e978034[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Nov 4 17:32:05 2020 +0100

    add statistical analysis for group analysis in diversity

[33mcommit 10845311559da2bd055c86508e4490895c8142f5[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Nov 4 16:13:58 2020 +0100

    Clarified some naming.
    count_starting_reads is a property!
    Finished draft of main downsize_pair function.

[33mcommit a91d584b9d444bb4dbd2770430589a7b54161174[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Nov 4 16:08:32 2020 +0100

    Renamed downsize_pair to ..._uncompressed.

[33mcommit d648867d1866f1d8e1c0a9ab5d19b99a324e24e2[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Nov 4 15:50:44 2020 +0100

    Removed extraneous IF statement based on filetype.
    This is no performed in read_downsize.py

[33mcommit 13cf8b553dc464541c655ad2e22344956f84ce65[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Nov 4 15:46:38 2020 +0100

    Added a property function to determine file type in the read_downsize.py
    module.
    
    Renamed the uncompressed downsizing function.
    
    Started new general downsizing main function that uses the file_type
    property to invoke the appropriate downstream function.

[33mcommit e0737e8c6488ab2d2d4d0adb8b47aa032e518fbe[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Nov 4 14:41:25 2020 +0100

    Corrected mistaken line downgrading Pandas.
    Restored to v1.0.1.

[33mcommit 55b76612772350cfee3d6f7b0d5b709379d9870e[m
Merge: d31a58b e81f93d
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Nov 4 10:57:30 2020 +0100

    Merge branch 'master' into 48-beta-diversity

[33mcommit e81f93dc6a316814ced8ac9d051379a1a05bc533[m
Merge: 94f1d1d d7408a1
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Wed Nov 4 10:16:49 2020 +0100

    Merge pull request #53 from motleystate/51-more-statistical-tests
    
    Statistical tests for group comparison

[33mcommit d7408a18b326df07b4935a36107e50727f0f51f4[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Nov 4 10:04:51 2020 +0100

    light renaming

[33mcommit 550f6462cc112ae03c8ea9c711cdb1a144e11f5a[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Nov 3 17:34:27 2020 +0100

    test conditions rather than using globals

[33mcommit 51dda40a7e0c9f711afc08c37438c260bf6e4df5[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Nov 3 14:29:04 2020 +0100

    delete df_split + remove eval

[33mcommit 794e3145d818290c3c114893f9d5f3e349a97c79[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Oct 30 16:45:49 2020 +0100

    Downsizing instance *kargs now passed from a dictionary, as a preliminary
    step towards running downsizing with Pool.

[33mcommit 4cd4b694aa92a6e2ebbbc51100cc7e19c5004d08[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Oct 30 15:52:02 2020 +0100

    little corrections

[33mcommit db4de9b96fe9599d221add559e8358d8fc0730a9[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Oct 30 15:24:37 2020 +0100

    Bug fixed where only starting read number, instead of all information, was
    passed to downsizing module.

[33mcommit cc242632564cbe2e358916dac4deddaf191d8e52[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Oct 30 13:49:43 2020 +0100

    Added the package *filetype* to requirements.txt for proper
     functioning of downsizing modules.

[33mcommit 7512615e1a6469c67fbb6b91cc14c9d4696e7bc0[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Oct 30 13:43:50 2020 +0100

    New count_starting_reads function, and changes to standard and gzip downsizing
     to use this generalized means of counting initial reads.

[33mcommit e361c3dfc2953027227e7a78e125a7b6fbdb1890[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Oct 30 11:29:16 2020 +0100

    output series or dataframe

[33mcommit 527b108fb7be322e4c0896530af99718bbe939e2[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Oct 29 11:29:59 2020 +0100

    ttest + chi2 contingency

[33mcommit 166213158efb1cf55ebc72825e408995b65a8192[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Oct 28 17:19:13 2020 +0100

    Moved file check from each downsize function to init.

[33mcommit 3408e9308ad04aea45b3962baaf04d59509c5f7f[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Oct 28 16:50:52 2020 +0100

    Moved output directory check to init.

[33mcommit c43b5e53880f5e8197cc160703ae1c866a0cbd5b[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Oct 28 15:57:09 2020 +0100

    Added some annotation to the MP code.

[33mcommit 3f6f7381fe40880f93ccba181304d98531ee25a5[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Oct 28 15:39:49 2020 +0100

    WIP! Not yet tested!
    
    -Added 'Pool' multiprocessing to down_dir_pair function, removing the existing
    single-process loop in read_info function.
    
    -Reformatted the resulting *list of dictionaries* to the expected dictionary.
    
    -Added option 'processes' to Class instantiation with some safeguards against
    running more processes than CPUs.

[33mcommit 94f1d1d5139879b0addd3092ac06ff53d4a3694a[m
Merge: 44fbd53 e63b49a
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Oct 26 17:01:22 2020 +0100

    Merge branch 'master' of github.com:motleystate/moonstone

[33mcommit 44fbd535e3de89d4e991d7d2350db3bedba3f845[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Oct 26 16:57:17 2020 +0100

    Adds gzip support to pair-ends reads downsizing and closes #33.
    Completes work to downsize all files in a directory and closes #34.

[33mcommit bf68f6f262ff362b3f757c9f7f0f42845213045a[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Oct 26 10:38:04 2020 +0100

    warnings empty group

[33mcommit d31a58bb382f0af6cbeb047c953bb4fe9fec3c1b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Oct 23 19:19:28 2020 +0200

    add scatter in graph

[33mcommit 3fe327e67422e28c48cef637efe498d2be9a592a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Oct 23 16:04:42 2020 +0200

    start pcoa for beta diversity

[33mcommit 4c860972e947c75b4d1a8d924d6e871fdc55ae9f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Oct 23 14:34:47 2020 +0200

    allow retrieving of grouped_df when analysing with metadata

[33mcommit 08fc8d025b5f74a84a39d0028e11eab55a6538b0[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Oct 23 12:15:41 2020 +0200

    flake8 corrections

[33mcommit ac4d9a01fac4b27973447af2a7492ee4e5cb0974[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Oct 23 12:14:06 2020 +0200

    allow customization of points to display for violin and box

[33mcommit da739aa28db4f5f4246b4aac675114d4e6c4e325[m
Merge: 86f7e36 e63b49a
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Oct 23 12:08:17 2020 +0200

    Merge branch 'master' into 45-mann-whitney-u-test

[33mcommit 86f7e36ab86cc45f7e019fa4e32ac9d07f3e0597[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Oct 23 12:05:27 2020 +0200

    mann - whitney-u

[33mcommit 535a9deb9a2d17d8469a5b3fc1263b1718a5100c[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 22 18:29:46 2020 +0200

    update documentation

[33mcommit e4191cc54bd55d53793600d0500efed589b903cc[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 22 18:23:06 2020 +0200

    refacto and add visualization for BetaDiversity

[33mcommit 09386721da82c6b755512771f221b91b3d10092b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 22 16:05:13 2020 +0200

    add BrayCurtis for beta diversity

[33mcommit e63b49ae0a2183fbba6909199f26ded9ab6d6d94[m
Merge: d2f8681 245b8c4
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Thu Oct 22 15:47:34 2020 +0200

    Merge pull request #44 from motleystate/43-more-filtering
    
    More filtering methods

[33mcommit 245b8c4db31aa4620ec95fe60b45ac0fe80e6fd5[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Oct 22 15:25:00 2020 +0200

    more tests

[33mcommit 47c3010cc10578fb564ba1855fa089268bcc84cf[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Oct 22 15:06:04 2020 +0200

    bug correction

[33mcommit e13dbeff1b53062bb15d3ad633fb96494db33b4e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 22 14:49:40 2020 +0200

    add test failing because of confusion in percentage naming

[33mcommit b54af6b86f8df047a2caa8bb8fd030267fc48146[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 22 14:37:51 2020 +0200

    start beta

[33mcommit aa2e046b71563ef7bad7b1b8bfbb57e6a770afb0[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Oct 22 14:20:17 2020 +0200

    update doc

[33mcommit 05a5bc3984b1f9f840ccc0d7138be747eb2022e7[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Oct 22 13:51:36 2020 +0200

    corrections

[33mcommit 079023a84d09cde525226d0b77c00d8044046db4[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Oct 22 10:26:59 2020 +0200

    df_split (2/2)

[33mcommit d2f8681bc6c1aa5cc259b226da210f5958454c69[m
Merge: 7084f14 f5d308b
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Wed Oct 21 17:46:24 2020 +0200

    Merge pull request #46 from motleystate/feature-box-plot-alpha-diversity
    
    Feature box plot alpha diversity

[33mcommit f5d308bd24a1fe73029eefbf365ec770b7707683[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Oct 21 17:26:22 2020 +0200

    update docstring for alpha indexes

[33mcommit fe47fa4c9e87a71707f5a509f497624962cff505[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Oct 21 17:15:03 2020 +0200

    handle log_scale in plotting_options

[33mcommit 418c64ba635f0b86a229b8a45e242b6ce38a6070[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Oct 20 18:38:21 2020 +0200

    add (log) to yaxis legend when applying log scale

[33mcommit dea73b94e26fb89edd21bc605814f745dfb9c511[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Oct 20 17:59:39 2020 +0200

    fix labels for GroupBaseGraph

[33mcommit 2729b376873cd062e6c3c4095f178c134eaea2f2[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 20 17:34:16 2020 +0200

    prerequisite df_split (1/?)

[33mcommit 5bab60bed3d982b50ecba0f5366cba3850e22b5a[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 20 16:33:26 2020 +0200

    filtering by number of different values

[33mcommit b0348e3162729d3217dfa7394f80ae1604bac574[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Oct 20 16:24:07 2020 +0200

    add option to force colors for groups

[33mcommit 920f975900fe2cafc2f24362441febe1c49028b6[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 20 14:55:03 2020 +0200

    filtering by percentage of nan

[33mcommit 16223e775e983d3f4832963893d7ca9aefa35ea1[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Oct 19 18:29:43 2020 +0200

    fix typo

[33mcommit cd4c814af9e1f33c68e1ad328d6521c3aa654cc4[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Oct 19 18:25:13 2020 +0200

    add boxplot visualization

[33mcommit 7084f141621bbeeb910c96c15563050059383368[m[33m ([m[1;33mtag: v0.4.0[m[33m)[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Oct 9 15:12:52 2020 +0200

    bump version for release

[33mcommit d2a0bac8e4a82acca4b25f1c76e4ba9d3f370f50[m
Merge: 58d1d7f 79302ea
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Fri Oct 9 12:15:44 2020 +0200

    Merge pull request #42 from motleystate/37-visualize-alpha-index-group-violin
    
    Visualize alpha index using violin plots

[33mcommit 79302ead11b136207056e3f0f88fcd14d3aa868a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 8 16:19:24 2020 +0200

    add docstrings for documentation

[33mcommit ae4f91e0535520f709468f96fe17b82ee03c0929[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 8 16:03:13 2020 +0200

    refactor way to plot groups with violin

[33mcommit 895c9ce24f3af332660ceec8669767630f291936[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 8 15:08:42 2020 +0200

    fix tests and flake8

[33mcommit 72bfd7f2d3173b75771e707df6bcf5d57723e52d[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 8 14:37:29 2020 +0200

    start violin visualization for alpha indexes

[33mcommit 9aa75ea483421902c9131a747431414b311755c9[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Oct 8 11:11:53 2020 +0200

    add way of defining index column while parsing metadata

[33mcommit 58d1d7fabc0ccd2ba3646029e20d5b8f0520ad21[m
Merge: 2525bbe c21dc8f
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Thu Oct 8 10:32:34 2020 +0200

    Merge pull request #41 from motleystate/40-sum-counts-in-reindex_with_taxonomy
    
    BUG : summing counts reindex_with_taxonomy

[33mcommit c21dc8f18e33bb2f79ae64a122f0b208fb02beaa[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 6 17:25:45 2020 +0200

    update docstring 2

[33mcommit 50a02e7c0ca623e5d324a028eb56d91bbc41d129[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 6 17:09:43 2020 +0200

    update docstring

[33mcommit 466106e6c9e93f12d613af065e6ca36ebe34f7cc[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 6 17:02:27 2020 +0200

    bug correction

[33mcommit ecaec1eb1840b06b8ab315d809c3100ed7a55616[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Oct 6 15:01:45 2020 +0200

    bug summing counts reindex_with_taxonomy

[33mcommit 2525bbe49d66b04eb5178b30b879f1236c9ec78d[m
Merge: da1b0a3 7a7d9d4
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Fri Sep 25 18:11:21 2020 +0200

    Merge pull request #36 from motleystate/35-handle-virgo-taxonomy
    
    Transform counts df into taxonomy df

[33mcommit 7a7d9d491e853c65fc5dee82e305c1f237ceb487[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 25 16:21:04 2020 +0200

    change print to logger.info

[33mcommit e5924e8dbd61d70dbb7833bd3018abc9545e6e5a[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Sep 25 12:30:32 2020 +0200

    Added 'downsize_dir' module to the project under normalization.reads
    Will call Downsize pair for each paired-ends read set.
    Todo:
     -Implement output directory in reads/base and read_downsize
     -Handle actual downsizing of gzipped files

[33mcommit 98c5d415ae96d35ffa0c927bfe3d28187a2ad25f[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 25 12:21:18 2020 +0200

    get rid off taxonomy_file

[33mcommit 59db0d5db060a166a50d8de11e751f4e979cfbf5[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Sep 24 17:18:10 2020 +0200

    GenesToTaxonomy (2/2 for now)

[33mcommit e76d859199e7fc72f150fdf04d714be5c7c118b2[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Sep 23 17:50:20 2020 +0200

    WIP - GenesToTaxonomy

[33mcommit cc14fe05b4399672dbce71075776155c2c862193[m
Merge: 228ad3b da1b0a3
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Sep 23 16:27:27 2020 +0200

    Merge branch 'master' into 35-handle-virgo-taxonomy

[33mcommit 228ad3bde7d1fb243a6ee8430a7d8264d9168a3b[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Sep 23 16:26:20 2020 +0200

    (2/2) reorganization taxonomy

[33mcommit da1b0a301986ab70381716bce02330c52546c3b0[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 23 10:24:45 2020 +0200

    hotfix: change char to concat for TaxonomyRandomSelection

[33mcommit 3b72d03c4bfea36f719dbfd80167ebf96d0f4281[m
Merge: 76fb767 88a1486
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Tue Sep 22 15:48:23 2020 +0200

    Merge pull request #32 from motleystate/31-randomselection-multiindexed
    
    31 randomselection multiindexed

[33mcommit 9c6e2918a598c1352ce9c29a79d1b20c51006637[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Sep 21 18:42:24 2020 +0200

    WIP - reorganization taxonomy

[33mcommit 0f61f73ebc3d5545387ce9dfc6230126e590f943[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Sep 21 14:16:12 2020 +0200

    removed extraneous print statment.

[33mcommit 76fb767745a238ec68cfb385318bfa1b688a18d1[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Sep 21 14:07:11 2020 +0200

    Corrects bug where downsized reads were not properly written.

[33mcommit 3feb4759e759128b1896498f280891f9c73aa2a6[m[33m ([m[1;33mtag: v0.3.0[m[33m)[m
Merge: 633d572 f5f3d48
Author: skennedy <seanpkennedy@gmail.com>
Date:   Fri Sep 18 16:34:45 2020 +0200

    Merge pull request #30 from motleystate/27-random-reads
    
    27 random reads

[33mcommit 88a148631189be09779ed081a34841ddfe477810[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Sep 18 15:15:33 2020 +0200

    add test to handle float in df

[33mcommit ec599fc0b2489af1f3320f23a67b0c1482ddeaf0[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Sep 18 11:56:19 2020 +0200

    update documentation

[33mcommit c978cb8e840a1321b1ed45a11430bb67a1ed0fb4[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Sep 18 11:55:28 2020 +0200

    add class to randomly select in multiindex df

[33mcommit f5f3d48ecbbd79f1d1431a9c6778c448124cf7f4[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Sep 18 08:26:28 2020 +0200

    Refactored [renamed] variables to more clearly indicate that they refered to file handles and not the actual data.
    Removed some redundant lines.

[33mcommit b154134e0d3427e19c6f0cf33a6e17ab592cd4b5[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Sep 18 07:48:52 2020 +0200

    Added missing closing of files for downsized reads.
    Changed print statments to logger.info

[33mcommit 914d65b23088fa9188729548e8dc8c2826002248[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Sep 17 17:31:58 2020 +0200

    fixed misplaced return in base module.
    Fixed variable assignment to accept tuple from downsize_pair function.

[33mcommit 633d572856cb10fe0ba142ebc75e4987cbb86f69[m
Merge: a850046 f2cb654
Author: skennedy <seanpkennedy@gmail.com>
Date:   Thu Sep 17 13:16:06 2020 +0200

    Merge pull request #28 from motleystate/26-randomly-select-counts
    
    26 randomly select counts

[33mcommit 00fda7da132e2d6081ac4138b36f6c1f5d51bd38[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Sep 17 13:03:29 2020 +0200

    Added a preset seed to randint.

[33mcommit 257be4a737ae694ee8cef9d4ab5f93953a36dd48[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Sep 17 12:54:53 2020 +0200

    Finished with paired-ends downsizing.
    Need to document usage.
    Still need to completed single read downsizing.

[33mcommit f2cb6546381789fdc45e658ed4cf4188ba001c42[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 17 12:39:26 2020 +0200

    add log for RandomSelection

[33mcommit a8500463b5a1c86a124b53e5440421c069e1b912[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 17 09:08:13 2020 +0200

    add numpy requirements to doc for scikit-bio

[33mcommit 0275633bdd705cf1604e321a51d458098468f966[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Sep 16 17:41:01 2020 +0200

    Added a __init__ file to the normalizations.reads directory.

[33mcommit 4320eaed6aaf0a98e2fd8c7d34a59cc04c3fbf65[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Sep 16 16:26:06 2020 +0200

    Finished draft of downsizing for paired end.
    Needs to be tested.

[33mcommit acbfdfc58398f2c9364489c7723a4de5024425c9[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 16 16:11:15 2020 +0200

    add documentation for normalization

[33mcommit ac1633eab7a39b7cd04b9ed346723cbb6a8cc0a8[m
Merge: a551351 0dab103
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 16 15:28:07 2020 +0200

    Merge branch 'master' into 26-randomly-select-counts

[33mcommit a55135131d555245a5ad43bd8e4930a1496b3f34[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 16 15:27:44 2020 +0200

    filter out samples below threshold

[33mcommit 0dab1039e0cd3257a7e9c77dde83ae2899ebb3f2[m
Merge: c157c71 87cf111
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Wed Sep 16 14:45:37 2020 +0200

    Merge pull request #29 from motleystate/25-alpha-diversity-Shannon-index
    
    Alpha diversity Shannon and Inverse Simpson indexes

[33mcommit 863ad92cd60eb45db0f98a475ba069ae4b6d2bee[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Wed Sep 16 11:40:53 2020 +0200

    Finished draft of downsizing base module.
    Renamed class and functions to 'downsizing'
    Decision made to seperate single and paired-ends reads.

[33mcommit 87cf11168e9e1fefaf9850902fb4721b76531727[m
Author: Agnes Baud <agnes.baud@pasteur.fr>
Date:   Wed Sep 16 10:47:03 2020 +0200

    update doc

[33mcommit 52e0dc6c81813b31df9ec617e66d9dc0036e6797[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 15 17:57:28 2020 +0200

    class SimpsonInverseIndex

[33mcommit 6326d7b8f544c7419dd8624b0b01f7a3b32da9b4[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 15 17:51:32 2020 +0200

    tests

[33mcommit aee704eaa8d643d7981cdaa38dcb300d741d2399[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 15 15:21:27 2020 +0200

    slight refactoring

[33mcommit 0ea0b79ce11a45e18b6d618660a557520ea2c06f[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 15 10:34:22 2020 +0200

    install numpy before

[33mcommit 621b6c6b4af038ed448062c0b2345df3e05de78f[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 15 10:02:57 2020 +0200

    update setup
    
    add scikit-bio to dependencies

[33mcommit 14d5982afee8280de6b62d0069c1a0c52ea4fb81[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Sep 14 18:04:20 2020 +0200

    class ShannonIndex

[33mcommit 4ca416e5388ca3530ea07bcb7b5e03ce4ef90abb[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Sep 14 17:23:12 2020 +0200

    minor changes after first use of plot module

[33mcommit d168d6de4826982159f1c2dead1165cf466b9887[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Sep 14 17:01:07 2020 +0200

    add class RandomSelection to normalization module

[33mcommit 9ec00630b58ab6f0cf21b24943c557280ce7a3f4[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Sep 14 16:02:33 2020 +0200

    Started work on base downsizing
    Added base module and pair-end downsizing module.

[33mcommit 1f4bc38ca2953a89245313a6d28049e753cbf1f1[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Sep 14 14:55:37 2020 +0200

    use BaseModule and BaseDf for BaseNormalization

[33mcommit c157c715e0295b5c1a3f46eb39abea1215cfc511[m[33m ([m[1;33mtag: 0.3.0[m[33m)[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 17:22:00 2020 +0200

    Changelogs and version

[33mcommit 22e29c3332e1e576521bc6b6bdd5e4461138705a[m
Merge: 75ee836 ffda325
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 17:01:39 2020 +0200

    Merge pull request #15 from motleystate/12-refactor-use-of-plotting-options
    
    Refactor use of `plotting_options`

[33mcommit ffda3251b44f35a8adda41cf3a8f4313f05b277d[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 16:59:17 2020 +0200

    analysis in index

[33mcommit a5216339093d76204b365b9faee93f8acf0caed5[m
Merge: 0880dd0 75ee836
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 16:56:58 2020 +0200

    Merge branch 'master' into 12-refactor-use-of-plotting-options ???

[33mcommit 75ee8364fef3ae1325c08f8f8973844fdbbce61c[m
Merge: 6664d9c 1afd8ab
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Fri Sep 11 16:48:51 2020 +0200

    Merge pull request #24 from motleystate/21-analysis-docs
    
    add DifferentialAnalysis to API doc

[33mcommit 0880dd0851ec54dd53d9d0a5a7fdc43e7dde775f[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 16:45:04 2020 +0200

    little doc update

[33mcommit 1afd8abbbe6bca5330ca5d75abe94e9dbf050691[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Sep 11 16:40:59 2020 +0200

    add DifferentialAnalysis to API doc

[33mcommit 5dc855a8313f306a0ac53e70f6b1343b9577253a[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 16:08:19 2020 +0200

    correction flake8 + 1 in heteregeneous_bins

[33mcommit a0abe43108c44275dd9f4d1d1cb2d809f6019507[m
Merge: 3e46ea5 6664d9c
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Sep 11 15:49:26 2020 +0200

    Merge branch 'master' into 12-refactor-use-of-plotting-options

[33mcommit 6664d9cbc2ac245a0e3259b65a9b2c2a631b008e[m
Merge: 33ddd4f 9b5f1d0
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Fri Sep 11 12:17:51 2020 +0200

    Merge pull request #17 from motleystate/13-BarGraph-refactoring
    
    BarGraph refactoring

[33mcommit 9b5f1d08cd610264ded124a87c012b5efb85ede8[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Sep 11 11:26:35 2020 +0200

    move binning for plot into utils.pandas

[33mcommit 6d9749c462aef16f66a542398ed69536e569f30e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Sep 11 10:45:21 2020 +0200

    fix distribution for plot sex in PlotMetadata

[33mcommit 3e46ea59f1635da4bae5813eaa1be450cd299c86[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Sep 9 17:46:38 2020 +0200

    more flexible plotting options

[33mcommit fec697f74bdffe6f1298140b9eafdf636714d9e7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 9 14:55:41 2020 +0200

    fix path in documentation

[33mcommit 2fd4342cd01a9b910629f54c15f9b626032f6098[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 9 14:28:29 2020 +0200

    refactor graph and plots

[33mcommit eb591f5a3b4d853dd63370eb596a4eb1bef2a9fb[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 9 14:01:16 2020 +0200

    get rid of tests of plotting options format and rename

[33mcommit caf3c5e39177f04e53d8eb0e5a2bf676aee448a2[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 9 12:24:26 2020 +0200

    fix tests for global plots prior to refactoring

[33mcommit 3366e442dc8764b158d69532f320e82791f5eb45[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 9 11:58:25 2020 +0200

    update plot templates

[33mcommit 64e439966a4bc331a386e9fc3852b9c18eaded33[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 8 17:50:41 2020 +0200

    log xaxis pass pytest

[33mcommit cdaa11c414f8f959d15fde4b7e9c1ef4f4a5c29c[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 8 17:33:10 2020 +0200

    log xaxis

[33mcommit 15f29c962c8e84de8bed42c735cb19da77cc06ae[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Sep 8 15:39:27 2020 +0200

    shapes in plotting_options

[33mcommit 19c75fb98e1eb695db07382e97aabe490893908b[m
Merge: 2eec553 33ddd4f
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Sep 8 10:59:03 2020 +0200

    Merge branch 'master' into 13-BarGraph-refactoring

[33mcommit 33ddd4f038653984678df161c60121aacbc71b4a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Sep 7 16:21:48 2020 +0200

    add __init__.py in plot module to fix #18

[33mcommit c9480628dacf782469abb8ceed7aec1151b04d1f[m
Merge: 09155f2 383826f
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Mon Sep 7 15:05:51 2020 +0200

    Merge pull request #19 from motleystate/14-available-classes
    
    14 available classes

[33mcommit 383826ff65a3a41626e978265f539fb79b3e703e[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Aug 31 17:44:02 2020 +0200

    Use of autosummary

[33mcommit 334ee1f9c05c69ee77f4547fb4a5afd3bc4af832[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Aug 31 13:16:21 2020 +0200

    plot available classes

[33mcommit 2eec553daf2acdc6cd46783e88e0c92d0d998792[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Aug 25 15:55:13 2020 +0200

    more tests for functions.py

[33mcommit 0f25e0270dff7b7e1b25eb85e6a41dcef1710e18[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Aug 24 18:09:04 2020 +0200

    independent BarGraph + some tests

[33mcommit 6934b2171ee9f13e89e025b37621802b04c63b20[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Aug 21 11:21:36 2020 +0200

    handling plotting options + show in plot_one_graph

[33mcommit 09155f2e1db077e7de7302a410e14db5344e4c4e[m[33m ([m[1;33mtag: v0.2.0[m[33m)[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Aug 20 16:48:17 2020 +0200

    update version for release

[33mcommit 88eda41fe6944b9a9ebd35888aca270d4ff426c8[m
Merge: 153d284 0fcccea
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Thu Aug 20 16:33:53 2020 +0200

    Merge pull request #9 from motleystate/3-plot-module
    
    plot module

[33mcommit 0fcccea8b727e4182ac008250ebb0b2b66edbbdb[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Aug 20 16:16:23 2020 +0200

    renaming to PlotCountsStats

[33mcommit ad2f9212535c3fbfedc211b211bb38c49f972b36[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Aug 20 16:03:21 2020 +0200

    update getting_started

[33mcommit 153d2847e816831c30cf9dcb9343f3d30369f248[m
Merge: cfeac82 f051085
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Thu Aug 20 15:35:41 2020 +0200

    Merge pull request #10 from motleystate/8-refactor-filtering-module
    
    Refactor filtering module

[33mcommit 1322a8ffd1bead9fea9dfe41c9b5bec349985a85[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Aug 20 14:19:36 2020 +0200

    abstractmethod decorator

[33mcommit f051085936ca1e021f0f2fa76886f17cd2839da2[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Aug 20 14:17:54 2020 +0200

    add validation for level for TaxonomyNamesFiltering

[33mcommit 07124202163381c78fc5e18617f8f2be8ea56b1a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Aug 20 14:07:26 2020 +0200

    add list of classes in documentation

[33mcommit 5976119129bec60610f0422223c5f8ecd33cfbcb[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Aug 20 10:21:05 2020 +0200

    update docs and changelogs

[33mcommit d75c72d265f250708c06b487ee29e1750603241f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Aug 19 19:19:09 2020 +0200

    remove old Filtering classes and update docs

[33mcommit 6328ae31ae225eb46846e621adfce9c11e425a71[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Aug 19 18:49:28 2020 +0200

    first corrections of comments

[33mcommit c870352a34f903a6b15b721ac5e48d267d3143eb[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Aug 19 17:37:08 2020 +0200

    fix typo and give explicit error message

[33mcommit 820c2cce17524bfac37598f7454b12c96d1b6288[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Aug 19 14:51:51 2020 +0200

    update plot

[33mcommit 96324098ceecbff5f1b671001ca074efbf135ea7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Aug 19 13:41:29 2020 +0200

    code gardening

[33mcommit 8ff7ab293471942dc067f6c4f7c684f29c31ef2d[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Aug 19 11:25:44 2020 +0200

    fixes for 3.6 ... hopefully

[33mcommit 67b331b8dc1bdc81f34878d013ef00ae5a5af2b0[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Aug 19 11:04:06 2020 +0200

    adapts typing type according to python version

[33mcommit 79c4292c489d184f94f0a05dd451f8fada176000[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Aug 18 20:02:11 2020 +0200

    docstring TaxonomyNamesFiltering

[33mcommit 0ae5d862d584fcabc5debfacd9dd839fa307ada9[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Aug 18 19:58:14 2020 +0200

    add TaxonomyNamesFiltering

[33mcommit d16dc9e969b3422da28a1389149fad36be5f411c[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Aug 18 18:29:12 2020 +0200

    more plot + checking plotting options (unfinished)

[33mcommit bdc54cada03e6217e83af62c08be5f1331150e9b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Aug 18 15:37:19 2020 +0200

    simplify exclusion of rows and columns with drop

[33mcommit a867036967edd7bc870450152e094d40a0d66323[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Aug 17 18:10:36 2020 +0200

    Add NamesFiltering that does not handle MultiIndex df

[33mcommit 1abe23ec0d1f621be92b1018746caceacd969723[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Aug 17 13:57:10 2020 +0200

    add way to either filter on column or rows for NoCountsFiltering

[33mcommit 562d3ed2a0e806bbfdd55ae7f81c0a929b85cc72[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Aug 17 13:36:00 2020 +0200

    change name and start extracting Filtering classes

[33mcommit 4edbf1ca129bfe4b21622851d0b9061adf1cd012[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Aug 14 18:32:03 2020 +0200

    update main bases and BaseFiltering

[33mcommit 3ab50971ae6180597e3e38e68dba310fddbd7dac[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Aug 14 16:58:19 2020 +0200

    using plotly instead of matplotlib

[33mcommit bd4775d43449eb4abf04e7d2d8dfaa0632f8f792[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Aug 13 18:02:47 2020 +0200

    add plotting_options

[33mcommit 0634f6ef355a88804c30999db7ee2f31f30a6518[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Thu Aug 13 08:19:20 2020 +0200

    organizational model

[33mcommit cfeac8298dce8f15d889ba9a582d8216e9be984b[m
Merge: c404ee4 bb6c206
Author: Agnès Baud <agnes.baud@pasteur.fr>
Date:   Wed Aug 12 11:48:06 2020 +0200

    Merge pull request #7 from motleystate/5-base-module-report-visualization
    
    base module for report and visualization

[33mcommit bb6c206dd5d84d3fec807ee21a38698492ae10c2[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Wed Aug 12 11:34:44 2020 +0200

    removal steps etc.

[33mcommit f2296a787741b75b6a86eb19a039ba9d777a2e1d[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Aug 12 10:43:20 2020 +0200

    small reformat of docs

[33mcommit d295612b29d5bd03b32650949d1b941a857e21ff[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Aug 11 16:05:46 2020 +0200

    resolves threshold to 0

[33mcommit 8ab778f8a12979a53f6d183149abb2439106890d[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Aug 11 15:32:22 2020 +0200

    doc for filtering_by_mean

[33mcommit c404ee4b38dc57d33f5713731d96ded074e7d087[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Aug 11 11:19:13 2020 +0200

    add pip-tools to requirements-dev.txt and add requirements.txt

[33mcommit 9de2246cb32a265f37a1c7832dcba295a0eb806c[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Aug 11 10:42:20 2020 +0200

    add failing test with threshold to zero for MeanFiltering

[33mcommit 3302a3cb033f6c96cb6180b4d0577b8bcdc27127[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Aug 10 17:35:16 2020 +0200

    renaming to mean_filtering.py + U base report_data

[33mcommit 7a0370740c8ba0a82d9e729e292e3f58def20b7e[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Aug 10 17:00:10 2020 +0200

    report_data + re-organisation (cf comments) etc.

[33mcommit 78725b0856b514a9182cedae7f104a1d1d0da69f[m
Merge: 5a36fad 4261410
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Aug 3 14:20:58 2020 +0200

    Merge branch 'master' into 5-base-module-report-visualization

[33mcommit 4261410a345f66132c32f27154e9d3e6dcf5007a[m
Merge: 9c551a7 1903cd3
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Mon Aug 3 09:28:23 2020 +0200

    Merge pull request #6 from motleystate/docs
    
    Start documentation with readthedocs

[33mcommit 5a36fad0d5b091637dd35d42c92b00147ea0231d[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Fri Jul 31 16:11:39 2020 +0200

    filteringbase + meanfiltering

[33mcommit 1903cd3deed11a31698792850f13270b9693baaf[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Jul 30 18:16:59 2020 +0200

    add doc about use of parsers

[33mcommit 2fafb612505fec1db979ea05f13acbf8a25949a3[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Jul 30 17:48:04 2020 +0200

    update documentation for parsers

[33mcommit 77e73d63fe82df460b54899b98f8e37357c37bdc[m
Author: Agnes Baud <agnes.baud@pasteur.fr>
Date:   Thu Jul 30 08:03:58 2020 +0200

    visualize for filtering module - 1st save

[33mcommit b2ac8e62d82912f20a1d98d06b0472e2d2af8e10[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jul 29 15:19:14 2020 +0200

    update documentation

[33mcommit da3e7a1b9572c4802acda9db0f96e25892747db4[m
Merge: e6dc10a 9c551a7
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jul 29 12:42:21 2020 +0200

    Merge branch 'master' into docs

[33mcommit e6dc10a49cff3ea8b7c1956a69a7519d0e668ce8[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jul 29 11:36:51 2020 +0200

    update changelogs in docs

[33mcommit 67ce96ab280304e1d5cd91c2a9c6969397e6b194[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jul 28 18:15:16 2020 +0200

    create module base in moonstone core

[33mcommit 8bf0761618c302c068da41407232cc051302edf9[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Jul 28 18:02:39 2020 +0200

    GetStarted page

[33mcommit 9c551a733ac36c99e65bc0f84d57b1f6a0edea7f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jul 28 17:28:38 2020 +0200

    add template for issues and pr

[33mcommit 0785439a6a0fc850c8a191f528b98825ce00c7ee[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Tue Jul 28 17:15:43 2020 +0200

    Create LICENSE

[33mcommit 5e153235e9e7f6c4df94ec66715d3a6016b791b9[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Tue Jul 28 17:14:37 2020 +0200

    Create CODE_OF_CONDUCT.md

[33mcommit 3c756087e8635a442266bab93ab2f8608caa53ae[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jul 28 17:04:00 2020 +0200

    add codecov badge

[33mcommit 21ed84a3a3db081317f69a578827ce560cf14e55[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jul 28 16:57:58 2020 +0200

    add coverage github action

[33mcommit 48bb2ecb6df61fb9538d44361f264de996f4e46e[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Tue Jul 28 13:35:59 2020 +0200

    add uninstall + list of dependencies

[33mcommit 88ec6c77be6b94c6b69bc95cf9698b03a42cb974[m
Author: Agnes BAUD <agnes.baud@pasteur.fr>
Date:   Mon Jul 27 18:07:12 2020 +0200

    unicode spe + start installation doc

[33mcommit 3ddd3c131b8b76ca24bb1f84c5d0562a3f6473fc[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 16:29:13 2020 +0200

    start structure for installation process docs

[33mcommit 2d6fc3d5640ded286e111822159fc62e9b35608a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 16:14:44 2020 +0200

    fix rtd issue

[33mcommit cce54f819fd31499428548deaf5ade0c8f2fa19f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 16:08:22 2020 +0200

    update index of documentation

[33mcommit 6872d18b30d1f1624362e9585b9cdf40a9d133ef[m
Merge: afe4662 547e1dc
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 16:05:49 2020 +0200

    Merge branch 'master' of github.com:motleystate/moonstone

[33mcommit afe466232b5f7cb7356e69af372eaac2a34e4381[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 16:05:39 2020 +0200

    start base for readthedocs documentation

[33mcommit 547e1dc0ed75b3b7c665b58dbb5066acfe7404d6[m[33m ([m[1;33mtag: v0.1.0[m[33m)[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Mon Jul 27 14:34:03 2020 +0200

    set up automatic publication on release

[33mcommit ab8e7b504512c81543357dc7d7c8bc1cc44f9f14[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 14:27:38 2020 +0200

    versioning and update README

[33mcommit 89f5010f5f9c431944da46563448a6ddeaf0343c[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 14:04:09 2020 +0200

    remove python 3.8.0 from CI github actions

[33mcommit bf396722fb1df5a19761d57b778858c681fb5587[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jul 27 13:56:57 2020 +0200

    remove python 3.5 from CI and rm gitlab ci

[33mcommit 212a17d14fd8e0bdff7e36ef7f07f9f5efd71592[m
Author: Kenzo-Hugo Hillion <hillion.kenzo@gmail.com>
Date:   Mon Jul 27 13:48:17 2020 +0200

    set up CI tests with github actions

[33mcommit 51164c9354dcc8b133048aede91085c9a0762d75[m
Merge: 83a2326 d4533b0
Author: Agnes  BAUD <agnes.baud@pasteur.fr>
Date:   Fri Jul 24 09:33:10 2020 +0200

    Merge branch 'abaud' into 'master'
    
    Add prefiltering module
    
    Closes #53
    
    See merge request metagenomics/data-analysis!25

[33mcommit d4533b026cd052bfa1cb0e5a695065356dddf8a2[m
Author: Agnes  BAUD <agnes.baud@pasteur.fr>
Date:   Fri Jul 24 09:33:10 2020 +0200

    Add prefiltering module

[33mcommit 83a2326b2e01bb22b634c07f74bea773c62d00eb[m
Merge: c262a4e 880c5eb
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Thu Jul 2 16:18:42 2020 +0200

    Merge branch '52-docs-and-readme' into 'master'
    
    update README
    
    Closes #52
    
    See merge request metagenomics/data-analysis!24

[33mcommit 880c5eb91fb1e4ac729ebd1ab92c169bacfd9305[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Jul 2 15:03:40 2020 +0200

    update README

[33mcommit c262a4ee79786501645d65df6b3ecdfdf5b21de5[m
Merge: e939854 78c1074
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jun 8 18:15:37 2020 +0200

    Merge branch 'iss51' into 'master'
    
    Issue 51
    
    Closes #51
    
    See merge request metagenomics/data-analysis!23

[33mcommit 78c1074ec0d55aa051ec8b9c8ece5c5821c64d0e[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Jun 8 17:38:53 2020 +0200

    one more test fix....fingers crossed!

[33mcommit a231b6d1b817523572c4bafc1cb6faaf9adb4d73[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Jun 8 17:33:33 2020 +0200

    attempt at fixing second test.

[33mcommit 21c9f426579eece6d960d3e5f7d6e4299d9a8ddf[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Jun 8 17:18:02 2020 +0200

    Try to correct failed tests. Scaling factors are now slightly different.
    Small text change.

[33mcommit fddbcce1e6463fdfab78b9fa1527d0b1f8aec222[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Jun 8 14:58:19 2020 +0200

    Too strict zero-filtering results in invalid scaling factors.
    Added progressive easing of zero_threshold parameter to computation of scaling factors.

[33mcommit 88e5bb28454a54a55b46e44ccd9bac4841198e57[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Jun 8 13:55:53 2020 +0200

    Fixed a bug that appeared under condition where zero filtering thresholds were low.
    This resulted in medians of zero, and invalid scalling factors.
    
    Removed replacement of NaN with Zero while computing scaling factors.
    Median now is based on numbers >0 only.

[33mcommit e939854f462db1f19456417214ceb1ff10f25228[m
Merge: 63b69b2 bbe89b3
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Mar 16 11:12:56 2020 +0100

    Merge branch 'merge_pop' into 'master'
    
    Merge pop
    
    See merge request metagenomics/data-analysis!22

[33mcommit bbe89b39247023bf6ba7c082958a387d5d12293e[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 15:57:34 2020 +0100

    In classify module:
    removed depricated merge function.
    fixed any broken merge calls in other functions (I hope!)

[33mcommit 6804f9757b579a05a2659d4c5fdcd9045b1127cc[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 15:42:27 2020 +0100

    finished test_df_merge. Test passes.

[33mcommit 79a93bf443f4ea75142566dfe0e1e40f0a6c6c71[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 15:19:01 2020 +0100

    started work on test_df_merge

[33mcommit 0bf34562c2e3bb91d9885ea5ff586ea1ffab64d6[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 14:24:39 2020 +0100

    added df_merge.py under /utils

[33mcommit 63b69b2563bf392922073dd229220659330caf9f[m
Merge: 767525b f3eff4d
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Mar 13 14:05:36 2020 +0100

    Merge branch 'ML_refactor' into 'master'
    
    ML_refactor_scaling
    
    See merge request metagenomics/data-analysis!21

[33mcommit f3eff4dba3dbff994fba04f3ca5274bd8fdd9c27[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Mar 13 14:01:38 2020 +0100

    updating names

[33mcommit c09b140a5495c336e61eaf304dea692be70f9180[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 13:21:53 2020 +0100

    removed unnecessary if __name__ == '__main__': in test_scaling_normalization.py

[33mcommit b6725601c769231275821b30353510f743014286[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 13:17:03 2020 +0100

    small change to correct calling @property in TestStandardScaler class

[33mcommit ead0062260bb902157bf7c28bbd569e72c3ddfe0[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Fri Mar 13 12:48:47 2020 +0100

    Corrected generic naming of scaling test issue.
    Changed expected type to NumPy array in scalTestStandardScaler

[33mcommit 9ba43a6d7dff6c1ccba3ce23f467bb552e6e29e8[m
Merge: 82cf2f4 767525b
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Mar 12 15:56:58 2020 +0100

    Merge branch 'master' into ML_refactor

[33mcommit 82cf2f4e3f66dfb95259ea8adfa357ecef19dc6b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Mar 12 15:56:54 2020 +0100

    fix flake8 and remove idea files

[33mcommit 3ad01ab5d13f06f0a40f0a71d6c28b3d422ad865[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Mar 12 11:56:18 2020 +0100

    start fixing flake8

[33mcommit d0375359093b35a9c4131c0e07804f75a9c5d865[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Mar 12 10:36:49 2020 +0100

    minor changes to classify.py

[33mcommit 18d97f7892f996112fb689d9d3404bda027e5560[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 10 15:57:32 2020 +0100

    Small tweaks to get the SVM module working.

[33mcommit dba725546415ba9504382cf5e86f36acb69e853b[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 10 15:36:31 2020 +0100

    changed SVM references to ML to reflact class name change

[33mcommit e11869cec872acf162cffeda66fd61fa963026fd[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 10 14:43:49 2020 +0100

    BaseScaler written
    StandardScaler first draft
    SVM in 'clasify' modified to use new scaler

[33mcommit 767525b42f744d46d69afd0b294b6e65a8dd3388[m
Merge: 3474df3 9d08ed9
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Wed Mar 4 08:18:37 2020 +0100

    Merge branch '44-refacto-qiime2-parser' into 'master'
    
    Refactor Qiime2 parser
    
    Closes #44
    
    See merge request metagenomics/data-analysis!20

[33mcommit 9d08ed9b774a6eb166b650218c2c19ae22c142a0[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Mar 3 14:44:11 2020 +0100

    refactor Qiime2 parser and tests

[33mcommit 11f21bbff15e25fb5b7f1e84324ff3fe79c5494a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Mar 3 14:40:11 2020 +0100

    add option to remove terms from taxa annotation

[33mcommit 82958f687c9ee5763fa0827ad1eae6e4f39a1274[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 3 14:26:05 2020 +0100

    Added NumPy v 1.18.1 to list of requirements.

[33mcommit 9915ed564ba60cecac93da3c1b6ebd8bfb0c18ad[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 3 14:25:18 2020 +0100

    Created new test module for scaler module

[33mcommit 6cf1c036076f29829be7dc88544ba27830bb640e[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 3 14:24:42 2020 +0100

    Created new ScalingNormalization module

[33mcommit bb6c6cee1c822e8fdf63274c492c948a8ce31acb[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Tue Mar 3 14:23:57 2020 +0100

    Simplified logging in module

[33mcommit 3474df378119078712278fa65a198289b55f17df[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Mar 3 11:55:38 2020 +0100

    rename test for geometric mean

[33mcommit 49180081bbf3e4fda43a1465a0a686cf049bb5a0[m
Merge: 13cdd65 c11f8b7
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Mar 3 11:53:55 2020 +0100

    Merge branch 'master' of gitlab.pasteur.fr:metagenomics/data-analysis

[33mcommit 13cdd65404a2824ef17b40366fe8619ef03e0a25[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Mar 3 11:53:44 2020 +0100

    add parameter to replace all 0 by 1 on input df for geometric mean normalization

[33mcommit c11f8b76bbbfeda98b4c000dca718865b83cd601[m
Merge: f24807a fc86a09
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Tue Mar 3 11:45:00 2020 +0100

    Merge branch '50-normalization-total-counts' into 'master'
    
    Total count normalization
    
    Closes #50
    
    See merge request metagenomics/data-analysis!19

[33mcommit fc86a09e5c8986d8a94408be038a8e60f0f81ffb[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Mar 3 10:59:40 2020 +0100

    add total counts normalization

[33mcommit f24807a6c7671fe8675b2eca3cef780708de4472[m
Merge: f755015 077e9b2
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Mon Mar 2 13:45:19 2020 +0100

    Merge branch '49-pandas-v-1' into 'master'
    
    Resolve "update to pandas v.1"
    
    Closes #49
    
    See merge request metagenomics/data-analysis!18

[33mcommit 077e9b28ea3922a37d923f31bb406347c43e4885[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 27 13:26:08 2020 +0100

    fix flake8

[33mcommit f89d868a6501f33217e0376e274397dd4a69dcfe[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 27 13:24:20 2020 +0100

    fix multi index path in pandas

[33mcommit 6930786868aa78bf13abc4b487bc76f4b17f3b9d[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 27 12:30:23 2020 +0100

    fix test in filtering on unexisting index

[33mcommit f7550159d8121d019cb71aed1ebf8bd642623b02[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 27 11:42:57 2020 +0100

    remove notebook files

[33mcommit 3360a89bb1d32dc316e33e6763ce28c198213a61[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 20 18:11:14 2020 +0100

    hotfix: deal with genus None species content

[33mcommit cbf07b031d6049f2496e7a0a6bc6a56d9e4b9a3f[m
Merge: 7b35677 17c3a08
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 20 18:06:33 2020 +0100

    Merge branch '46-fix-species-name' into 'master'
    
    Merge genus and species when specified
    
    Closes #46
    
    See merge request metagenomics/data-analysis!17

[33mcommit 17c3a08ac53d2373ca1b649f6225a37d22d20a09[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Feb 20 18:04:29 2020 +0100

    merge genus and species when specified

[33mcommit 7b3567704db67d8edc26bc5acca89e0a266636c5[m
Merge: 6ad0a16 763ad6f
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Thu Feb 20 10:31:59 2020 +0100

    Merge branch '43-update-normalization' into 'master'
    
    Change the way normalization level is handled
    
    Closes #43
    
    See merge request metagenomics/data-analysis!16

[33mcommit 763ad6fdeb9a43fa64b6b11be2d65b89bb36a07e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Feb 19 15:02:39 2020 +0100

    cleaning code for normalization

[33mcommit 741e55258aaa59422d85bfb14f349756d27b7959[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Feb 19 11:53:47 2020 +0100

    fix flake8

[33mcommit ea10994c4341deffc141aa86a240526929292cd4[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Feb 19 11:50:59 2020 +0100

    no normalization level by default

[33mcommit 4ef85522472274f28d38331c943076462fec7523[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Feb 19 11:26:52 2020 +0100

    change name of normalization

[33mcommit 6ad0a162c51f3bb78f9c45d4044be14f33bb7ab2[m
Merge: ae59b07 da6b955
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Wed Feb 19 10:05:20 2020 +0100

    Merge branch '41-kraken2-merge-parser' into 'master'
    
    add sunbeamkraken2parser
    
    Closes #41 and #42
    
    See merge request metagenomics/data-analysis!15

[33mcommit da6b95575094b400519dc9d4c6f3e33c0e92388e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 15:52:59 2020 +0100

    finish refactor for kraken2 and metaphlan2 parsers

[33mcommit d66e884378a871aa762cf46e55555c0da39ab81f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 15:43:27 2020 +0100

    Start refactoring for #42

[33mcommit 5e47fbab2191dceb7348780325e1b06d663c0eee[m
Merge: 749d2b7 ae59b07
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 15:34:56 2020 +0100

    fix conflict

[33mcommit ae59b0744a7f778ddd34321ba099f8dcda5962d7[m
Merge: c74f4c2 d365adc
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 15:33:59 2020 +0100

    Merge branch '40-parser-metaphlan2-merged' into 'master'
    
    Parser for metaphlan2 merged table
    
    Closes #40
    
    See merge request metagenomics/data-analysis!14

[33mcommit d365adc6776343c71ab3693ad0c8fe1211ddb053[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 15:29:25 2020 +0100

    add test for t__ example

[33mcommit 749d2b741f3ad6faecde7f5af1821de8a68c1bc0[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 15:18:55 2020 +0100

    add sunbeamkraken2parser with tests

[33mcommit 8b72428abc0fc58535cd1ced4ac5fb39f7f6de4f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 14:08:03 2020 +0100

    fix tests

[33mcommit 86e9e5c939b435cf65d2c4c967076f8d54b5fbd7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 14:03:34 2020 +0100

    deal with t__ from metaphlan2

[33mcommit 135797031adbbaa587c21f1cbb7015918b56598b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 13:25:38 2020 +0100

    fix flake8

[33mcommit 48849780a52d4455b97c491891cfd4c3c8ee2fa8[m
Merge: 220374e c74f4c2
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 13:23:41 2020 +0100

    Merge branch 'master' into 40-parser-metaphlan2-merged

[33mcommit 220374e2e937812466e9d7564befe5ede8926e82[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 13:22:05 2020 +0100

    add dataframe with multiindex

[33mcommit b36b048253176fdc17f7d7a5783a3c92139cb6b8[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 13:09:43 2020 +0100

    add method to split taxa column to different columns

[33mcommit c74f4c282a24e0fe0fd4c515e7fccc11ab39c97a[m
Merge: 36e4263 2504bed
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Tue Feb 18 10:41:10 2020 +0100

    Merge branch 'rfecv-dev' into 'master'
    
    Redevelopment of 'Recursive Feature Elimination' (RFE)
    
    See merge request metagenomics/data-analysis!9

[33mcommit 2504bedcc12f6c583bde1f87b85df47ba900d51a[m
Merge: aff413b 0502f17
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 09:48:30 2020 +0100

    Merge branch 'rfecv-dev' of gitlab.pasteur.fr:metagenomics/data-analysis into rfecv-dev

[33mcommit aff413b2180623e5151d2aa463d46cfd7c717df7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 09:48:21 2020 +0100

    force versions of pandas to 0.25.3 and scikit-learn to 0.21.3

[33mcommit bccb643417d6baeccdb72bd43ab6541ad7553a13[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Feb 18 09:32:11 2020 +0100

    start metaphlan2 parser

[33mcommit 0502f17e466b41e8f31a97ce8dd7a832a01fe0b7[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Mon Feb 17 16:03:52 2020 +0100

    module classify: figures now saved. logging more extensive. ready to merge.

[33mcommit 36e4263ef8985512bdba1ddda03a30e437e56fc2[m
Merge: ae5487c 18f272f
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Wed Jan 29 18:06:35 2020 +0100

    Merge branch 'cleaning-normalization' into 'master'
    
    Reorganize normalization module and tests
    
    See merge request metagenomics/data-analysis!13

[33mcommit 18f272f6fe7d01b8c99670ae18cdfe35bc71a26e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jan 29 17:38:25 2020 +0100

    Add .df to be current state of dataframe during normalization

[33mcommit 5b5c615b66d57b76370e83cbb9e9dcb400d60d25[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jan 29 16:40:52 2020 +0100

    Reorganize modules to correspond to architecture

[33mcommit a81f83d0b183ca1b062996790c8945fac4a5ed5a[m
Merge: cff9bc6 ae5487c
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jan 29 10:06:41 2020 +0100

    Merge branch 'master' into cleaning-normalization

[33mcommit ae5487c7a1884a2fa57d3e9ddec74936ec3f1a01[m
Merge: 0c59209 c65cc71
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Wed Jan 29 09:10:17 2020 +0100

    Merge branch 'tidying-up-parser' into 'master'
    
    Tidying up parsers
    
    See merge request metagenomics/data-analysis!12

[33mcommit cff9bc622818652da06375a2eb74a3838d8fce75[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 28 17:07:44 2020 +0100

    reorganize normalization module

[33mcommit c65cc7156404cbb0c470f4677c0be798f6288ece[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 28 15:00:06 2020 +0100

    make the standardize dataframe the one by default for Qiime2 parser

[33mcommit afd564aab210b6bf3cdd4f6e32e91eccc33e61d4[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 28 13:33:10 2020 +0100

    update organization of code for parsers and transform

[33mcommit 168d5d9e909b1007f3ee8a3dbfd9f438305c5cb4[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 28 10:25:14 2020 +0100

    add DataFrameCleaner to transform imports

[33mcommit a7ab2a812a666dd82b9b15cd7f4b9f949a1c22f7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 28 10:24:12 2020 +0100

    Move transform operations to parsers module

[33mcommit 246af5a0ab6a4afabfe2aeca4451bf09e6ba6425[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Jan 24 18:02:45 2020 +0100

    reorganize parsers

[33mcommit b9fe402dd61b11fd172961b2b16f98cd53b32e8c[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Jan 23 17:02:45 2020 +0100

    update README and add contributors

[33mcommit efa8fdcd98129ca2d832fa866f87ba1400137aa5[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Jan 16 15:31:54 2020 +0100

    Added additional logging for merging dataframes.
    Save SVC coeficients graph.
    Display RFECV plot.

[33mcommit 0c592094141e66b1ded50675e61240afed8a9046[m
Merge: d5d8c01 b82984d
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Wed Jan 15 10:04:50 2020 +0100

    Merge branch '27-yaml-config-cleaning' into 'master'
    
    Apply cleaning while parsing based on YAML config file
    
    Closes #27
    
    See merge request metagenomics/data-analysis!11

[33mcommit b82984de3fd9781b8f9f5a77cf8a483d2c127c13[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 14 11:50:52 2020 +0100

    explain in more details the example given in README

[33mcommit 7c5dd63fb4c0b306a1a5b3a701440f7d64c2411d[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Jan 14 11:46:47 2020 +0100

    add few words on the parsing based on YAML

[33mcommit 3b5e93631e37cb4a652d3720fc90234586863036[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jan 13 17:35:01 2020 +0100

    add end to end test for basic usage from YAML config file

[33mcommit 0b5eb81a4310cbfc2ce93f9b63918ff4224217b7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jan 13 15:28:34 2020 +0100

    add rename step in parsing test

[33mcommit ac722f7aaa62b76451ad3e3186633ff34315c55b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jan 13 15:21:14 2020 +0100

    add renaming method and change historization

[33mcommit 2dd21e00de371672e4b87589f40169e5a6950dbd[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jan 13 14:55:48 2020 +0100

    add cleaning steps in the Metadata parser

[33mcommit 076ee9c6a88f886f362b59214c6508204b00575d[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jan 13 14:18:34 2020 +0100

    refactor parser and update depending code

[33mcommit 24103727c9d4f9489c14df7491e9eb76e369d205[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Jan 13 12:46:55 2020 +0100

    update history and moved test data

[33mcommit 41c7b4d7af6d4ee8e88e39fbe72f1c10b5a36b6b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Jan 10 17:47:17 2020 +0100

    add string cleaning

[33mcommit a2e4631c433aee533062ed006bf621b2ae3d772b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Jan 10 15:53:42 2020 +0100

    create base for transformations
    
    fix flake8

[33mcommit d5d8c01d5cb2e25dd207de6aa10284084596dc9a[m
Merge: 96247ed f9a6256
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Jan 9 10:27:45 2020 +0100

    Merge branch 'improve-metadata-statistics' into 'master'
    
    add value repartitions to stats number
    
    See merge request metagenomics/data-analysis!10

[33mcommit f9a6256aa037b472870ed07c6be6b6f563d8c17e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Jan 8 17:05:59 2020 +0100

    add value repartitions to stats number
    
    fix broken tests
    
    remove comments

[33mcommit 96247ed053977c11e45b0d247ef6aeda79ade138[m
Merge: 602ebb5 3e3c5c1
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Dec 19 15:23:31 2019 +0100

    Merge branch 'mergingdf' into 'master'
    
    Merging dataframes and differential analysis
    
    Closes #25
    
    See merge request metagenomics/data-analysis!8

[33mcommit 3e3c5c18ff697ad560f36131fbf8da6fee296474[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 18 17:39:15 2019 +0100

    changing subtle naming of a variable

[33mcommit fd83d49f77c8410324d767fdbf1dea640ef5f431[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 18 15:18:20 2019 +0100

    adding function to del only zeros rows and test

[33mcommit 0716199cd6daa9f4c52714fb5b29e995401cc847[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 18 14:58:28 2019 +0100

    Brief explanation  how  module works conceptually

[33mcommit 54f27aa456fda27a9ca5ac6d5b9c9c2eb4e90b21[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 18 14:51:10 2019 +0100

    building unique dichotomic test & multiple test

[33mcommit 5b3712e1eb8ca6ee600a537ff237a6c23cac57fa[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Dec 18 11:14:09 2019 +0100

    add logger

[33mcommit 4f878d446afcf4065be8aebd984340114ed12e6f[m
Merge: dfed155 602ebb5
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Dec 18 11:10:27 2019 +0100

    Merge branch 'master' into mergingdf

[33mcommit dfed155e4954a531d46e6ee1b71504b89dd7819a[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 18 11:11:39 2019 +0100

    changing read to read_counts

[33mcommit 8caf83505eeeeaddd2c1dbcd12c4439a69cfb19b[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Dec 17 16:47:20 2019 +0100

    adding the missing blank line

[33mcommit f1a510885119dd03b8904829961d402ac3339145[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Dec 17 16:45:08 2019 +0100

    Improving normalization module to keep multiindex

[33mcommit 6d5236bc28ca7130504582acb4c24f73b33f765d[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Dec 17 15:12:53 2019 +0100

    removing rows with only zeros

[33mcommit 946378befa42764e6b23f222b96e71796b597f0c[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Dec 17 10:48:19 2019 +0100

    correcting the class name

[33mcommit 65c7bded6165d7dd8c89ac97b66d7090ffe38bed[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Dec 16 17:36:55 2019 +0100

    finish tests for differentialanalysis

[33mcommit 2d1301ea095f7b83c7acaf6f2ee2c7342f74367b[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Dec 16 17:08:47 2019 +0100

    add another packages requirement in setuptools

[33mcommit a7e51c0311e9649cb2ddb18fb8b84a57daa5a872[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Dec 16 16:17:07 2019 +0100

    changing test names of test file

[33mcommit fe83ff5e6af7240825e3c14200c4369fd9eb3ce3[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Dec 16 16:13:20 2019 +0100

    restructuring diff analysis module, adapt old test

[33mcommit 602ebb5ae979e78ff48b8f20418cd3f3a0b6b0fa[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Dec 12 18:15:25 2019 +0100

    create class for statistics about df

[33mcommit fb45221ae925a7d9cd736824b219bc1012cbe61a[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Dec 12 17:18:48 2019 +0100

    quick fix metadata parser

[33mcommit 69db33892bee04b0350a8d2efebe995ce110f69e[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 17:04:03 2019 +0100

    deleting an old file version

[33mcommit a32ea0f964d648cd854fddc4104f29245a12591f[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 16:49:39 2019 +0100

    removing extra white line

[33mcommit 1398d5b0151a748d25084c1d43f30873f8c33a68[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 16:47:13 2019 +0100

    trying out again with good pytest raise excep

[33mcommit ceca5a93452c5d8556f3103d1a0ffd1aebe93e27[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 16:45:31 2019 +0100

    trying out thing  because it does not recognise change

[33mcommit a38f404964429bda488f4d4c9bbbc590eb6f6ab8[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 16:42:50 2019 +0100

    just to try

[33mcommit 6b92397e5b9586a77f963d8a3f82bd308a1ef42a[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 16:35:51 2019 +0100

    typo, full table position and raise exception test

[33mcommit 786fd11d744819a97c96758d88f304279293fad1[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Dec 11 13:25:49 2019 +0100

    fix import in test

[33mcommit 03f92a817d1ec7c204f2ed03afb290c52ec7e370[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 13:15:51 2019 +0100

    deleting the file

[33mcommit 2736397fa10853dfc386443e627123f57c1d75c0[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 13:12:44 2019 +0100

    changed folders and files

[33mcommit ff37059f5dec1f4a582e7d9dfc529aecc0818fd7[m
Merge: d238f15 407c952
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Dec 11 12:20:51 2019 +0100

    Merge branch 'master' into mergingdf

[33mcommit d238f15ca1e4b151bd76d035cc9eb3f473f06dd9[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 11:59:05 2019 +0100

    changing  module names and raising  error

[33mcommit 407c9526f550cf44ee7b716825f7d364bd34e27a[m
Merge: 43b9f21 694f3eb
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Dec 11 11:27:50 2019 +0100

    Merge branch 'filtering' into 'master'
    
    Filtering
    
    See merge request metagenomics/data-analysis!7

[33mcommit 694f3eb3c2e0eca49156c1c6e451b34f8b713808[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 11:21:04 2019 +0100

    changing module name and directory location

[33mcommit 7b3f38a15959d8b95638d405c04a3026151b5bd3[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 11:16:11 2019 +0100

    changing corrections in attributes

[33mcommit ce45b0605e20b587a65c411bcdf2262814379560[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 10:11:29 2019 +0100

    getting rid of the extra white space

[33mcommit b4aabe76982b946a00e70959ad0bec5a09672b88[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 10:09:52 2019 +0100

    adding explanation for self.step

[33mcommit cc3afb799f9d80ad79127aa08f71edbc65b3fc0f[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 10:04:35 2019 +0100

    removing white space

[33mcommit cd2f6ad43c2b4702d7f1f10d4eece42b1f2e47a3[m
Merge: 8e4c056 462dbf9
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 10:02:13 2019 +0100

    fixing some of the corrections

[33mcommit 8e4c056d8d75182f9ff6187c6e5c86ba38d7713c[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 10:00:53 2019 +0100

    changing part of the comments

[33mcommit 8ed601b5aa8884590a584278ca778ceca1989db4[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 09:45:48 2019 +0100

    getting rid of the extra white spaces at end line

[33mcommit ca4fa933c5c7447fab992f3c20cc8f29786f53e4[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 09:41:55 2019 +0100

    checking formating options

[33mcommit 190985dc180617696f5a1d3df9632ebccbf3ec7a[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 11 09:31:03 2019 +0100

    correcting  formating mistakes

[33mcommit 462dbf9a681ae64929b7817920a82057220dca5b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Dec 10 14:44:33 2019 +0100

    add test for filtering method without relevant info

[33mcommit 0a8374f8ccaca0cda6d19f6400c616b51c71bdbe[m
Author: Mariela F <mariela.furstenheim@pasteur.fr>
Date:   Tue Dec 10 13:19:31 2019 +0100

    tests for differential analysis

[33mcommit dcb2fa3f8a8fd6261899c81270635977ddac94c1[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Dec 9 15:01:24 2019 +0100

    differential analysis  and concat module (related)

[33mcommit 3c8b1dfbb3068539a06fbb65dc834581f3843254[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Dec 9 11:53:51 2019 +0100

    creating the concat class and tests

[33mcommit 22ee20a76076949c25e1ccf8c9bf20273f6b8254[m
Merge: 988e62c 43b9f21
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Dec 6 16:56:52 2019 +0100

    Merge branch 'master' into filtering

[33mcommit 988e62c81ac6901376a23897db381b6ef8bab009[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Dec 6 16:41:01 2019 +0100

    filter row & column based on metadat & option input

[33mcommit 43b9f2118300ce8c4e7a8dc0485aeaebd833d4b5[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Dec 6 10:39:15 2019 +0100

    fix flake8

[33mcommit 9dcc540145070dd726e7963798181904a27de5ca[m
Merge: 7605b1b 11f6158
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Dec 6 10:36:12 2019 +0100

    Merge branch 'df_standardisation' into 'master'
    
    Df standardisation
    
    Closes #19
    
    See merge request metagenomics/data-analysis!6

[33mcommit 11f615831752cd6de7900534e371783d87ee4550[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Dec 6 10:36:38 2019 +0100

    adding skipping multiple rows option

[33mcommit 91c42eb4615c60e7a800b7412e0fae2e91e35cf2[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Dec 4 14:38:23 2019 +0100

    modified filling missing taxa values with format

[33mcommit 1013c74c703872ec6131c8658a068f8cf049c74b[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 29 15:34:39 2019 +0100

    including nheritance

[33mcommit cb6b0d75379f17ddfa709c8d9dc69fc1d5dd994b[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 29 13:14:15 2019 +0100

    changing test location and importing series

[33mcommit 7605b1b3d07903e46be8fb2a5a4f5526ed6384c3[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 13:08:34 2019 +0100

    Saving SVM components to a file.
    Removed this print statement.
    
    Started work on DataAnalysis Type.

[33mcommit a322461f58cfdfadf5509d88528999c3c3f50c0f[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 29 12:29:47 2019 +0100

    reformated import test

[33mcommit dbbfdeba82f158e8b680a4be36c03eeaf2a9f2c2[m
Merge: 3393370 a39c050
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 11:19:06 2019 +0100

    Merge branch 'logger_integration' into 'master'
    
    Logger integration
    
    See merge request metagenomics/data-analysis!4

[33mcommit a39c05050e13921a930de764dce320b45a6b8856[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Fri Nov 29 11:06:10 2019 +0100

    fix flake8

[33mcommit 955f8e8315ed9c668008806eb9a689318d8ebfb3[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 10:58:03 2019 +0100

    Logger code added to the remaining modules.

[33mcommit 4604985f79c909ab639d8356fbdcb10a692f2087[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 10:49:51 2019 +0100

    Added logging to handleCounts.

[33mcommit b7bf759fa409e37e12393783b51d56c7cc0e691e[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 10:43:38 2019 +0100

    Generalized logging of Classes by fetching __class__.__name__
    
    Added logging for handleMetadata.py

[33mcommit 355ea9f3cdd85a57a806b5a8830ffbaf61d2ecfe[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 29 10:34:39 2019 +0100

    shortened the length to pass test

[33mcommit e3164ddac57d6b45653a1264f9c8013992f7a8ec[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 29 10:30:17 2019 +0100

    Suggestions changed and # in column name fixed

[33mcommit 9fa301fa898f931ecae11256f3238ed6ee4475bd[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 10:17:10 2019 +0100

    Logging now present in all functions in clasify.py
    
    Also saving PDF of ROC analysis.
    This required passing output directory to this module.

[33mcommit 1402b4cdd5982a18b6b8245cc3d41ee437aca610[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Fri Nov 29 09:03:16 2019 +0100

    Quick change to save random forest plot instead of opening a display window.

[33mcommit 44fb4c37385fab42418de0ae2587481f4266f287[m
Author: skennedy8 <seanpkennedy@gmail.com>
Date:   Thu Nov 28 18:55:29 2019 +0100

    Added logging structure for randomforrest.
    
    Need to handle report file and some user feedback.

[33mcommit 226e996251e19fe9e0689efb7af35ce7a6fb06c9[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Thu Nov 28 16:02:47 2019 +0100

    Added logging text to randomForest.
    
    Still need to remove print statments.

[33mcommit 3a2d9bfc64fdc3996743fdd6ebc4682858db9d7a[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Thu Nov 28 15:55:27 2019 +0100

    Generalized logging format with __name__ convention.
    Applied to 'filtering.py'

[33mcommit c89875499192774c543be64758691f4b5b0c747a[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Thu Nov 28 15:17:00 2019 +0100

    Silenced the matplotlib logging spam!

[33mcommit f1329b077a35c67c76b689514153200e2ad35336[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Thu Nov 28 15:10:00 2019 +0100

    Added the basics of the logger 'moonstone_app' in main.
    A file to record all logging levels.
    A stream for ERROR and above.

[33mcommit 903de862e82f9c3aceaf7cf4c7773c6a059e2635[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Thu Nov 28 14:09:07 2019 +0100

    WIP: need to remove # from OTU ID

[33mcommit e5dfe79dfb7b28069a1795d87f17cc04d30147cc[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Thu Nov 28 14:03:15 2019 +0100

    Corrected previous correction to target the FIRST column: index_col=0
    Renamed the index column of the metafile as 'sample'
    Small modifs in main to get output directory to the openmeta method, and save the imported metadata

[33mcommit 868beb3fc973c4abf54bd0b4de64cecd48b1e927[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Thu Nov 28 13:35:25 2019 +0100

    df_standardisation change into qiime & moved csv

[33mcommit 45d7432c2107a673402edf819634cd38ff573991[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Thu Nov 28 12:54:40 2019 +0100

    Fixed a bug whereby the index column of the metadata files had to be called 'samples.' If there was no such column, or if another so-named column exists, bad things would happen.
    
    This was changed so that it is simply the first column that is used as the index when importing.

[33mcommit eb8da2416d3981a20dd344430bab515c5cfa6be1[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Thu Nov 28 09:45:55 2019 +0100

    creating the qiime module

[33mcommit bd3cee44c9e463d7f407f1d1cf0bf6eb920a2fab[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 27 17:24:14 2019 +0100

    deleting the print statement

[33mcommit 33933709399bd3467846740d1a85f5420c10d149[m
Merge: 09c65fa 0dd559f
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Nov 27 16:43:59 2019 +0100

    Merge branch '9-Values_to_log' into 'master'
    
    Creating logger and non-zero threshold option
    
    Closes #18
    
    See merge request metagenomics/data-analysis!5

[33mcommit 0dd559f9454a5d398148ab8b5beaa377746737e7[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 27 12:34:58 2019 +0100

    removed blank line at the end of the file

[33mcommit 417292f52a5f77263ab286ae38965b23e875b340[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 27 12:32:31 2019 +0100

    small flake8 error solved

[33mcommit aac92d5ae7876103d896657e6ca58f94f1eb35a6[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 27 12:28:31 2019 +0100

    corrected test for removed zero df

[33mcommit 0956176dc328a1e29909c8627686ed2862ae3e80[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 27 12:20:24 2019 +0100

    adding the appropiate removed zero log function

[33mcommit 7b3647c940ac9f67b8ca7707ff2fa13e5453a0b3[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 26 14:58:12 2019 +0100

    Fixing a bug of the data(from str to int)

[33mcommit 2d9fcdc8c7321643db831dfda35c55f68979c942[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 26 14:40:16 2019 +0100

    changing complete filepaths to relative paths

[33mcommit 8db7be3ec29abeb23e85a0413220b209e2d4bac3[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 26 14:32:45 2019 +0100

    Functions commented with examples

[33mcommit 61f9b1dee33fc6eb52ff3c9881fd4db4ad2a498a[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 26 12:40:46 2019 +0100

    Df_standardisation completed and tests too

[33mcommit 68e57224ecfe7a7624e006fa0ffefbb019551985[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Nov 25 15:33:51 2019 +0100

    WIP:  df_standardisation

[33mcommit e9bef2b70c4f48d65b86294f11b55ae839ff1d16[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 22 09:36:57 2019 +0100

    adding test for zero_df and changed int in thresh

[33mcommit e8ff62718cf55e6e26301509ce13d49992793832[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Thu Nov 21 12:22:24 2019 +0100

    Including the limminting cases in normalization

[33mcommit 1b518485b79a99c25644ec7c444eb1aafb619864[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 21 10:11:02 2019 +0100

    fix test to have the expected behaviour

[33mcommit f9f6db3a763f6924371b4426f3022e183be5bf49[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 20 10:09:52 2019 +0100

    fixed test for the threshold of zeros in a row

[33mcommit 925cf8debd4e0d5f87210d6513c43e060b15cbd8[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 19 15:22:15 2019 +0100

    WIP:  accounting for a number of zeros threshold.

[33mcommit 09c65fad599d406ec1fc39bc7c4ba73bf47f698c[m
Merge: b97b6f2 da0d357
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Mon Nov 18 12:11:26 2019 +0100

    Merge branch '9-Values_to_log' into 'master'
    
    Resolve "Df normalization"
    
    Closes #1, #9, #10, #11, #12, #13, and #14
    
    See merge request metagenomics/data-analysis!2

[33mcommit da0d3579d4edae1160bf9165656cb20ef568b120[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Nov 18 12:10:05 2019 +0100

    changed file names and added init files in modules

[33mcommit 8a4c3b6c20e16c4075c4c7535117b1da27bc4783[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Nov 18 11:43:30 2019 +0100

    changed names of classes

[33mcommit c9fce072a14aa9a7e9771db927b837160c92f3cd[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Nov 18 11:15:58 2019 +0100

    changing the name of scaling factors

[33mcommit 8aaaee42f062cab853351a247f76afd1853491f2[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Mon Nov 18 10:36:56 2019 +0100

    applymap for logs and base exponents

[33mcommit 9fe5f1280dbb803ed4f23fc2155b8c1b1f6cec61[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Fri Nov 15 14:47:59 2019 +0100

    Missing how to change the base exponent in CSF

[33mcommit 8204f879e68559d04c7f6b20b118475385821441[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Thu Nov 14 12:00:15 2019 +0100

    WIP: changing dependance between methods.

[33mcommit 7afa281bf87c2845d973e3d5a5f2c2265a121b75[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Wed Nov 13 10:46:25 2019 +0100

    finish the tests and removed the extra properties

[33mcommit b97b6f2481fcbbdbd607592c43bf61108e94d53e[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Tue Nov 12 16:02:04 2019 +0100

    Starting adding logging to replace 'print' statements.
    
    Modified filtering.py only prior to creating a marge request.

[33mcommit 2b90977487d1e47db2c2949e07f874458374669c[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 12 15:18:49 2019 +0100

    Added most of the tests and including some properties

[33mcommit 3ccec7dc2bd9cf500a143fac192c38c7d846f74b[m
Merge: 1dd9c72 c1bea86
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Tue Nov 12 15:13:06 2019 +0100

    Merge branch '15-refactor-moonstone-script' into 'master'
    
    refactoring main script of moonstone
    
    Closes #15
    
    See merge request metagenomics/data-analysis!3

[33mcommit c1bea8607d22bbf85d7edcde0dcc4baf0a745c7f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Nov 12 14:53:31 2019 +0100

    Test and finish refactoring for main script

[33mcommit 1dcdf566de1acbc8920a64a3030efee21a709dac[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Nov 12 13:03:48 2019 +0100

    refactoring main script of moonstone

[33mcommit 875a5f36714e9b2c0ba8da3ed88ec4a58e2b3d90[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 12 12:58:04 2019 +0100

    adding two properties and renaming scaling_df func

[33mcommit 645afd7b7e7b0b640614208003410df847db00f6[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 12 12:47:35 2019 +0100

    Erased importing pandas at the beggining

[33mcommit ee4dac7e034a152db0dfa18be03675ec957ec731[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 12 12:36:48 2019 +0100

    Added the functions left to obtain the scaled df

[33mcommit cd9b713e29074b8f059a1e9e520d0872b6e289ed[m
Author: Mariela <mariela.furstenheim@pasteur.fr>
Date:   Tue Nov 12 12:22:32 2019 +0100

    Remove zeros and Log values of the df

[33mcommit 1dd9c7262e58ce9ce13a94898c3de3d75af6548d[m
Merge: 2e3e3fb 39ab4a2
Author: Sean  KENNEDY <sean.kennedy@pasteur.fr>
Date:   Tue Nov 12 11:41:51 2019 +0100

    Merge branch 'importGitHUb' into 'master'
    
    Import skennedy code from github
    
    Closes #3 and #7
    
    See merge request metagenomics/data-analysis!1

[33mcommit 39ab4a21da70532f827def6246ac1c5f6b64db22[m
Author: skennedy <sean.kennedy@pasteur.fr>
Date:   Tue Nov 12 11:17:11 2019 +0100

    Final changes before merge to Gitlab.
    
    -Streamlined if/else for directory creation.
    -Added a default value to output directory argument.

[33mcommit 0d436f5d9aa6a0eee11d8e6d415d347d1c5ac295[m
Author: Sean KENNEDY <sean.kennedy@pasteur.fr>
Date:   Tue Nov 12 10:10:12 2019 +0100

    Directed output files to output directory:
    
    -starting_variables.csv (stats.py 32-39)
    -filtered_variables + filteredCountFile : both use stats.py. (also main 91-93)
    -metadata_withKClusters (clustering.py 86-88)
    -rf_AllFeatures (randomForest.py 122-123)
    
    Now passing output directory argument to Classes and the output filename to methods.

[33mcommit 2478d07e5be301fba3a31c442f8b53b4419a664c[m
Author: Sean KENNEDY <sean.kennedy@pasteur.fr>
Date:   Tue Nov 12 09:14:31 2019 +0100

    Added an dedicated output directory as a required option.
    
    In main:
    50-56: creates directory is it does not exist
    37-47: adds a safeguard against overwriting the existing directory. This warning can be suppressed with an option.

[33mcommit cc01c5a2d356c1a194886a2a78236f86bd577dea[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 7 16:51:02 2019 +0100

    make moonstone executable from command line

[33mcommit 13416070a106aa3372fcc2e5754769bfeaccb228[m
Author: Sean KENNEDY <sean.kennedy@pasteur.fr>
Date:   Thu Nov 7 15:53:09 2019 +0100

    Added all Python file for Moonstone.
    From GitHub.

[33mcommit 2e3e3fba54036fd33124ce6cb08b657766b13f45[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 7 14:20:25 2019 +0100

    update README with moonstone name

[33mcommit 315ae76ff7c69f97a5f7f9e03a505b397a0bdf09[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 7 14:19:12 2019 +0100

    remote notebook from repo

[33mcommit 30541e479c5daa61f287d1ea860d1ec8abb8677f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Nov 7 14:14:38 2019 +0100

    Change name to moonstone

[33mcommit 8d89cf727f9ffd0942c34f80c2d444e9415d1ff6[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 26 12:08:23 2019 +0200

    Update README

[33mcommit 27031e34522502ed5b06ac676727a3ca3757ea9d[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 26 12:03:01 2019 +0200

    add notebook example

[33mcommit ecd539772b93809a48b7bddf580114262fed302b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 26 10:25:57 2019 +0200

    Allow metadata without header

[33mcommit 7902e7b686b5f7671d6920c04c644d772f26b452[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 26 10:03:33 2019 +0200

    refacto and build class for Series stats

[33mcommit 966e25c7dc941c16ab8c1f63dcd6e37f093cf58e[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Thu Sep 26 09:38:36 2019 +0200

    add method to return stats of columns from dataframe

[33mcommit ae959cb2b358f0633d3c0e6dafea3a5e0bc04c80[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 25 16:30:04 2019 +0200

    add methods to calculate stats on dataframe

[33mcommit df38aeb3218b43deaefe391ed2a5a8c1b2ed4b7f[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 25 16:02:35 2019 +0200

    fix dataframe property for base

[33mcommit d5192f85232d88c0bc8cabbf9a44962c1a3d8f1b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 25 11:50:51 2019 +0200

    add converter from pandas dtype to python type

[33mcommit af75f71d0b12e8e874482428e98f9529ff7921c8[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 25 11:44:23 2019 +0200

    add function to print multi col elements

[33mcommit 08d3a73c14a61f60794d482214adda49ce6365f2[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 25 10:53:55 2019 +0200

    remove python3.6 from CI try to fix coverage
    
    debug CI for files
    
    cat .coverage file
    
    add __init__ to tests folder
    
    Restore python3.6 on CI

[33mcommit 7e38e90724f01f291586308af85da47d554bc4e2[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Sep 24 17:03:09 2019 +0200

    add coveragerc to skip tests and setup.py

[33mcommit 1a8229812979113299b8124e60568ce343ec48ce[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Sep 24 16:59:31 2019 +0200

    add badges

[33mcommit ffa6cd7d141bdd558605c442cb1cb94e8b1285bf[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Sep 24 16:57:12 2019 +0200

    add metadata parser

[33mcommit 40167204a2fdd333db4222da35a9bbe05cb27f8b[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Sep 24 16:31:52 2019 +0200

    fix ci

[33mcommit 535fa7a7aff9237367c6a31482fe04bbb2e049c7[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Tue Sep 24 16:30:58 2019 +0200

    add test for picrust2 parser

[33mcommit d9560200512391f494543894832927c714d9a6db[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 4 17:37:06 2019 +0200

    add parser for picrust2 results

[33mcommit 77fc3631c866cd002aeab9d88d046b1e568dd9df[m
Author: Kenzo-Hugo Hillion <kenzo-hugo.hillion1@pasteur.fr>
Date:   Wed Sep 4 15:30:19 2019 +0200

    Initial commit
