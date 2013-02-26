"""
The intial version of this code was adapted from a recipe by Zoran Isailovski
(published under PSF License).

http://code.activestate.com/recipes/413268-dependency-injection-the-python-way/
"""

import logging


class FeatureBroker(object):
    def __init__(self, allowReplace=False):
        self.providers = {}
        self.allowReplace = allowReplace

    def Provide(self, feature, provider, *args, **kwargs):
        if not self.allowReplace:
            assert feature not in self.providers, \
                "Duplicate feature: %r" % feature

        if callable(provider):
            call = lambda: provider(*args, **kwargs)
        else:
            call = lambda: provider
        self.providers[feature] = call

    def __getitem__(self, feature):
        try:
            provider = self.providers[feature]
        except KeyError:
            raise AttributeError("Unknown feature named %r" % feature)
        return provider()

    def __setitem__(self, feature, component):
        self.Provide(feature, component)


features = FeatureBroker()


def register(feature, component):
    """register `component` as providing `feature`"""
    features[feature] = component
    return component

## Representation of Required Features and Feature Assertions

# Some basic assertions to test the suitability of injected features

def NoAssertion():
    def test(obj):
        return True
    return test


def IsInstanceOf(*classes):
    def test(obj):
        return isinstance(obj, classes)
    return test


def HasAttributes(*attributes):
    def test(obj):
        for attr_name in attributes:
            try:
                getattr(obj, attr_name)
            except AttributeError:
                return False
        return True
    return test


def HasMethods(*methods):
    def test(obj):
        for each in methods:
            try:
                attr = getattr(obj, each)
            except AttributeError:
                return False
            if not callable(attr):
                return False
        return True
    return test

# An attribute descriptor to "declare" required features

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
    def __init__(self, feature, assertion=NoAssertion(),
                 alt_name="_alternative_"):
        self.feature = feature
        self.alt_name = alt_name
        self.assertion = assertion
        self.result = None

    def __get__(self, callee, T):
        self.result = self.Request(callee)
        return self.result  # <-- will request the feature upon first call

    def __set__(self, instance, value):
        '''Rename the feature'''
        if isinstance(value, str):
            logging.info("changed %s to %s" % (self.feature, value))
            setattr(instance, self.alt_name + self.feature, value)
        else:
            raise TypeError("can't change the feature name to non-string type")

    def __getattr__(self, name):
        assert name == 'result', \
            "Unexpected attribute request other then 'result'"
        return self.result

    def Request(self, callee):
        fet_name = self.feature
        if hasattr(callee, self.alt_name + self.feature):
            fet_name = getattr(callee, self.alt_name + self.feature)
        obj = features[fet_name]

        try:
            obj.register_handler(callee)
        except AttributeError:
            pass

        isComponentCorrect = self.assertion(obj)
        assert isComponentCorrect, \
            "The value %r of %r does not match the specified criteria" \
            % (obj, self.feature)
        return obj

class OptionalFeature(RequiredFeature):
    def Request(self, callee):
        fet_name = self.feature
        if hasattr(callee, self.alt_name + self.feature):
            fet_name = getattr(callee, self.alt_name + self.feature)

        try:    
            features[fet_name]
        except AttributeError:
            return None

        return super(OptionalFeature, self).Request(callee)

class Component(object):
    "Symbolic base class for components"
    def __init__(self):
        self.observers = []

    @staticmethod
    def _rm_duplicate_deps(deps):
        for i, d in enumerate(deps):
            if d in deps[i + 1:]:
                del deps[i]
        return deps

    def get_dependencies(self):
        deps = [o.get_dependencies() for o in self.observers]
        deps = sum(deps, self.observers)
        deps = Component._rm_duplicate_deps(deps)
        return deps

    def register_handler(self, handler):
        if handler not in self.observers:
            self.observers.append(handler)

    def unregister_handler(self, handler):
        if handler in self.observers:
            self.observers.remove(handler)

    def notify_observers(self):
        for dep in self.get_dependencies():
            dep._update()

    def _update(self):
        pass

    def update(self):
        self._update()
        self.notify_observers()


class dictproperty(object):
    """implements collection properties with dictionary-like access.
    Copied and modified from a recipe by Ed Swierk
    published under PSF license
    <http://code.activestate.com/recipes/440514-dictproperty-properties-for-dictionary-attributes/>`_
    """

    class _proxy(object):
        def __init__(self, obj, fget, fset, fdel):
            self._obj = obj
            self._fget = fget
            self._fset = fset
            self._fdel = fdel

        def __getitem__(self, key):
            if self._fget is None:
                raise TypeError("can't read item")
            return self._fget(self._obj, key)

        def __setitem__(self, key, value):
            if self._fset is None:
                raise TypeError("can't set item")
            self._fset(self._obj, key, value)

        def __delitem__(self, key):
            if self._fdel is None:
                raise TypeError("can't delete item")
            self._fdel(self._obj, key)

    def __init__(self, fget=None, fset=None, fdel=None, doc=None):
        self._fget = fget
        self._fset = fset
        self._fdel = fdel
        self.__doc__ = doc

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return self._proxy(obj, self._fget, self._fset, self._fdel)
