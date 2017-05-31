# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazo@correo.nucleares.unam.mx
r"""A script to obtain the percentage of passed doctests."""

from os import system
# from fast.config import fast_path

system("python tests.py > doctest_output.txt")

f = file("doctest_output.txt", "r")
lines = f.readlines()
f.close()

items_without_tests = []
tests_total = []
items_total = []
tests_failed = []
tests_passed = []
for line in lines:
    # We check for lines reporting how many items had no tests
    if "items had no tests" in line:
        n = line[:line.find("items had no tests")]
        items_without_tests += [int(n)]

    # We check for lines reporting how many tests there were in all items.
    if "tests" in line and "items." in line:
        ntests = line[:line.find("tests in")]
        nitems = line[line.find("tests in")+8:line.find("items.")]

        tests_total += [int(ntests)]
        items_total += [int(nitems)]

    # We check for lines reporting how many tests passed and how many failed.
    if "passed and" in line and "failed." in line:
        npassed = line[:line.find("passed and")]
        nfailed = line[line.find("passed and")+10:line.find("failed.")]

        tests_passed += [int(npassed)]
        tests_failed += [int(nfailed)]


coverage = (1.0 - float(sum(items_without_tests))/sum(items_total))*100
success = float(sum(tests_passed))/(sum(tests_total))*100

print "The coverage of doctests is", coverage, "%."
print "The sucess of the doctests is", success, "%."
