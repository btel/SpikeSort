from SpikeBeans import base
from nose.tools import ok_,raises
from nose import with_setup

def setup():
    "set up test fixtures"

def teardown():
    "tear down test fixtures"
    base.features = base.FeatureBroker()

class Dummy(base.Component):
    con = base.RequiredFeature('Data', base.HasAttributes('data'))

    def get_data(self):
        return self.con.data
    
class DummyDataProvider(base.Component):
    data = "some data"

class NixProvider(base.Component):
    pass

@with_setup(setup, teardown)
def test_dependency_resolution():
    base.features.Provide('Data', DummyDataProvider)
    comp = Dummy()
    ok_(comp.get_data()=='some data')

    
@raises(AssertionError)
@with_setup(setup, teardown)
def test_missing_attribute():
    base.features.Provide('Data', NixProvider)
    comp = Dummy()
    data = comp.get_data()
    
@raises(AttributeError)
@with_setup(setup, teardown)
def test_missing_dependency():
    comp = Dummy()
    data = comp.get_data()
    