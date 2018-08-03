
import sys, os
import shutil

from setuptools import setup, Command, find_packages

try:
    import yaml
except ImportError as e:
    print("Setup requires package 'pyyaml'")
    sys.exit(1)


# +-----------------------------------------------------------------------------
# | CONSTRUCT PARAMETERS FOR setuptools
# +-----------------------------------------------------------------------------

def trunc_lines(s):
    parts = s.strip().split('\n')
    parts = [p for p in parts if len(p) > 0]
    return ''.join(parts)

def build_keyword_dictionary(prefs):
    keywords = {}

    for key in [
        'name','license','url','download_url','package',
        'package_dir','platforms','description','install_requires',
        'long_description','package_data','include_package_data','scripts'
        ]:
        if key in prefs:
            keywords[key] = prefs[key]

    keywords['author'] = \
        ', '.join(prefs['authors'][:-1]) + ' and ' + \
        prefs['authors'][-1]

    keywords['author_email'] = \
        ', '.join(prefs['emails'])

    keywords["package_dir"] = \
        {package: '/'.join(package.split('.')) for package in prefs['packages']}

    keywords['long_description'] = \
        trunc_lines(keywords['long_description'])

    output = ""
    first_tab = 40
    second_tab = 60
    for key in sorted(keywords.keys()):
        value = keywords[key]
        output += key.rjust(first_tab) + str(value).rjust(second_tab) + ""

    return keywords


with open('setup.yml','r') as f:
    yaml_string = ''.join(f.readlines())
    preferences = yaml.load(yaml_string)


setup_args = build_keyword_dictionary(preferences)

setup (**setup_args)


