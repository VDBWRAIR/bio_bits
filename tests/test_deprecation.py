from . import unittest, THIS
import mock

@mock.patch('bio_bits.version.get_release', mock.MagicMock(return_value=1.0))
class TestProjectRename(unittest.TestCase):
    def test_can_import_bio_pieces(self):
        import bio_pieces
        r = bio_pieces.version.get_release()
        self.assertEqual(1.0, r)

    def test_can_import_sub_module(self):
        from bio_pieces import version
        r = version.get_release()
        self.assertEqual(1.0, r)
