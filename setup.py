from setuptools import setup

if __name__ == "__main__":
    setup(
        entry_points = {
            'console_scripts': ['rnachromprocessing=RnaChromProcessing.main:main',
                                'detect-strand=RnaChromProcessing.detect_strand:main',
                                'infer-xrna=RnaChromProcessing.infer_xrna:main'],
        }
    )
