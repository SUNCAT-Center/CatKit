import os
import unittest
import shutil
from catkit.hub.postgresql import CathubPostgreSQL
from catkit.hub import db2server, make_folders_template, folder2db

path = os.path.abspath(os.path.join(os.path.dirname(__file__)))

class UploadTestCase(unittest.TestCase):
    def setUp(self):
        """Set up a temp file."""
        if not os.path.exists('temp'):
            os.makedirs('temp')

        """Ensure postgres database is empty"""
        db = CathubPostgreSQL(user='postgres')
        con = db._connect()
        db._initialize(con)
        db.truncate_schema()

    def tearDown(self):
        """Clear temporary files."""
        shutil.rmtree('temp')

    def test_make_folders(self):
        template_data = {
            'title': 'Fancy title',
            'authors': ['Doe, John', 'Einstein, Albert'],
            'journal': 'JACS',
            'volume': '1',
            'number': '1',
            'pages': '23-42',
            'year': '2017',
            'publisher': 'ACS',
            'doi': '10.NNNN/....',
            'DFT_code': 'Quantum Espresso',
            'DFT_functional': 'BEEF-vdW',
            'reactions': [
                {'reactants': ['2.0H2Ogas', '-1.5H2gas', 'star'],
                 'products': ['OOHstar@top']},
                {'reactants': ['CCH3star@bridge'],
                 'products': ['Cstar@hollow', 'CH3star@ontop']},
                {'reactants': ['CH4gas', '-0.5H2gas', 'star'],
                 'products': ['CH3star@ontop']}
            ],
            'bulk_compositions': ['Pt'],
            'crystal_structures': ['fcc', 'hcp'],
            'facets': ['111']
        }

        title = template_data['title']
        authors = template_data['authors']
        journal = template_data['journal']
        volume = template_data['volume']
        number = template_data['number']
        pages = template_data['pages']
        year = template_data['year']
        publisher = template_data['publisher']
        doi = template_data['doi']
        dft_code = template_data['DFT_code']
        dft_functional = template_data['DFT_functional']
        reactions = template_data['reactions']
        crystal_structures = template_data['crystal_structures']
        bulk_compositions = template_data['bulk_compositions']
        facets = template_data['facets']

        make_folders_template.main(
            title=title,
            authors=authors,
            journal=journal,
            volume=volume,
            number=number,
            pages=pages,
            year=year,
            publisher=publisher,
            doi=doi,
            DFT_code=dft_code,
            DFT_functional=dft_functional,
            reactions=reactions,
            custom_base='temp',
            bulk_compositions=bulk_compositions,
            crystal_structures=crystal_structures,
            facets=facets,
        )

    def test1_read_folders(self):
        folder2db.main('{path}/aayush/'.format(path=path))

    def test2_upload(self):
        db2server.main('{path}/aayush/MontoyaThe2015.db'.format(path=path),
                       user='postgres')
        if os.path.exists('{path}/aayush/MontoyaThe2015.db'.format(path=path)):
            os.remove('{path}/aayush/MontoyaThe2015.db'.format(path=path))


if __name__ == '__main__':
    unittest.main()
