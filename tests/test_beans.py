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
    
    def __init__(self):
        self.data = False

    def get_data(self):
        return self.con.data
    
    def on_data_change(self):
        self.data = True
    
class DummyDataProvider(base.Component):
    data = base.DataAttribute("some data")

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

@with_setup(setup, teardown)
def test_on_change():
    base.features.Provide('Data', DummyDataProvider())
    comp = Dummy()
    comp.get_data()
    base.features['Data'].data = "new data"
    ok_(comp.data)
