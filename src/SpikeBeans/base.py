import os

######################################################################
## 
## Feature Broker
##  
######################################################################

class FeatureBroker:
    def __init__(self, allowReplace=False):
        self.providers = {}
        self.allowReplace = allowReplace
    def Provide(self, feature, provider, *args, **kwargs):
        if not self.allowReplace:
            assert not self.providers.has_key(feature), "Duplicate feature: %r" % feature
        if callable(provider):
            def call(): return provider(*args, **kwargs)
        else:
            def call(): return provider
        self.providers[feature] = call
    def __getitem__(self, feature):
        try:
            provider = self.providers[feature]
        except KeyError:
            raise AttributeError, "Unknown feature named %r" % feature
        return provider()


features = FeatureBroker()

######################################################################
## 
## Representation of Required Features and Feature Assertions
## 
######################################################################

#
# Some basic assertions to test the suitability of injected features
#

def NoAssertion(): 
    def test(obj): return True
    return test

def IsInstanceOf(*classes):
    def test(obj): return isinstance(obj, classes)
    return test

def HasAttributes(*attributes):
    def test(obj):
        for each in attributes:
            if not hasattr(obj, each): return False
        return True
    return test

def HasMethods(*methods):
    def test(obj):
        for each in methods:
            try:
                attr = getattr(obj, each)
            except AttributeError:
                return False
            if not callable(attr): return False
        return True
    return test

#
# An attribute descriptor to "declare" required features
#

class DataAttribute(object):
    """A data descriptor that sets and returns values
       normally and notifies on value changed.
    """

    def __init__(self, initval=None, name='var'):
        self.val = initval
        self.name = name

    def __get__(self, obj, objtype):
        return self.val

    def __set__(self, obj, val):
        self.val = val
        for handler in obj.observers:
            handler()


class RequiredFeature(object):
    def __init__(self, feature, assertion=NoAssertion()):
        self.feature = feature
        self.assertion = assertion
        self.result=None
    def __get__(self, obj, T):
        self.result = self.Request(obj)
        return self.result # <-- will request the feature upon first call
    def __getattr__(self, name):
        assert name == 'result', "Unexpected attribute request other then 'result'"
        return self.result
    def Request(self, callee):
        obj = features[self.feature]
        try:
            #handler = getattr(callee, ("on_%s_change" % self.feature).lower())
            handler = callee.update
            obj.register_handler(handler)
        except AttributeError:
            pass
            
        assert self.assertion(obj), \
                 "The value %r of %r does not match the specified criteria" \
                 % (obj, self.feature)
        return obj

class Component(object):
    "Symbolic base class for components"
    def __init__(self):
        self.observers = []
    def register_handler(self, handler):
        if handler not in self.observers:
            self.observers.append(handler)
    def unregister_handler(self, handler):
        if handler in self.observers:
            self.observers.remove(handler)
    def notify_observers(self):
        for handler in self.observers:
            handler()   
    def update(self):
        self.notify_observers()
            

######################################################################
## 
## DEMO
## 
######################################################################

# ---------------------------------------------------------------------------------
# Some python module defines a Bar component and states the dependencies
# We will assume that
# - Console denotes an object with a method WriteLine(string)
# - AppTitle denotes a string that represents the current application name
# - CurrentUser denotes a string that represents the current user name
#

