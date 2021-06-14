# https://github.com/shayfletcherz/NumericalAnalysisEx16.git

# EX16

import sympy
import sympy as sp
from datetime import datetime
from datetime import datetime

local_dt = datetime.now()
d = str(local_dt.day)
h = str(local_dt.hour)
m = str(local_dt.minute)

# Divide the given segment into small segments of 0.1
from sympy import lambdify

#Division of a segment into several segments according to parameter n
def split_segment2(a, b, n):
    distance = abs(b - a)
    counter = n
    segments_list = [a]

    for i in range(counter):
        a = a + (distance / n)
        segments_list.append(a);

    return segments_list

#Dividing a section into several sections with a length of 0.1
def split_segment(a, b):
    distance = abs(b - a)
    counter = int(distance / 0.1)
    segments_list = [a]

    for i in range(counter):
        a = a + 0.10
        segments_list.append(a);

    return segments_list

#Finds a sign replacement in a selected section
def find_sign_change(segments_list, function):
    sign_change_list = []

    for i in range(len(segments_list) - 1):
        if function(segments_list[i]) * function(segments_list[i + 1]) < 0:
            sign_change_list.append((segments_list[i], segments_list[i + 1]))

    return sign_change_list


def newton_raphson(a, b, function, derivative, final_func, epsilon):
    xr = (b + a) / 2
    counter = 1

    while counter < 100:

        print("Iteration number:", counter, "| f(x):", function(xr), "| f'(x):", derivative(xr), "| XR:", xr)
        xr1 = xr - (function(xr) / derivative(xr))

        if abs(xr1 - xr) < epsilon:
            break

        xr = xr1
        counter += 1

    if counter >= 100:
        print("The equation does not converge")
        return

    if 0 <= abs(final_func(xr)) <= 0 + epsilon:
        print("The number of iteration is :", counter)
        print("The root is :", str(xr) + '00000' + d + h + m)

def secant_method(a, b, function, final_func, epsilon):
    xr = a
    xr_1 = b
    counter = 1

    while counter < 100:

        print("Iteration number:", counter, "| f(xr+1):", function(xr_1), "| Xr+1:", xr_1)
        xr1 = (xr_1 * function(xr) - xr * function(xr_1)) / (function(xr) - function(xr_1))

        if abs(xr1 - xr) < epsilon:
            break

        xr_1 = xr
        xr = xr1
        counter += 1

    if counter >= 100:
        print("The equation does not converge")
        return
    r = final_func(xr)

    if 0 <= abs(final_func(xr)) <= 0 + epsilon:
        print("The number of iteration is :", counter)
        print("The root is :", str(xr) + '00000' + d + h + m)

    return


def hight(a, b, n):
    return (b - a) / n

def simpson(function, segment, h):
    results = [function(segment[0])]

    print("Iteration number:", 0, "|x:", segment[0], "|f(x):" ,function(segment[0]))

    for i in range(1, len(segment) - 1):

        if i % 2 == 0:
            results.append(2 * function(segment[i]))
            print("Iteration number:", i, "|x:", segment[i], "|f(x):" ,function(segment[i]))
        else:
            results.append(4 * function(segment[i]))
            print("Iteration number:", i, "|x:", segment[i], "|f(x):" ,function(segment[i]))

    results.append(function(segment[len(segment) - 1]))
    print("Iteration number:", i+1, "|x:", segment[len(segment) - 1], "|f(x):" ,function(segment[len(segment) - 1]))

    return 1 / 3 * h * sum(results)


# The function returns the area in the range by romberg method
def rombergMethod(function, startPoint, endPoint, limit, epsilon):
    results = [[0 for i in range(limit + 1)] for j in range(limit + 1)]  # creation of matrix
    for k in range(0, limit):
        res = trapezMethod(function, startPoint, endPoint, 2 ** k)  # calculate trapez method
        results[k + 1][1] = res  # storing values
        print("R" + str(k + 1) + "," + str(1) + " = " + str(res))  # print
    for j in range(2, limit + 1):
        for k in range(2, limit + 1):
            results[k][j] = results[k][j - 1] + (
                        (1 / ((4 ** (j - 1)) - 1)) * (results[k][j - 1] - results[k - 1][j - 1]))
            print("R" + str(k) + "," + str(j) + " = " + str(results[k][j]))  # print
            if abs(results[k][j] - results[k - 1][j]) < epsilon:  # check if the difference is less then epsilon
                return results[k][j]
    return results[j - 1][k - 1]

# function returns The area in the range by trapez method
def trapezMethod(function, startPoint, endPoint, segments):
    x = sp.symbols('x')
    function = lambdify(x, function)
    h = (endPoint - startPoint) / segments
    sum = 0
    while startPoint < endPoint:  # run unti end point
        sum += 0.5 * ((startPoint + h) - startPoint) * (function(startPoint) + function(startPoint + h))
        startPoint += h
    return sum


def main_methode():
    x = sp.symbols('x')
    function = (x**2 * sp.exp(-x**2+5*x-3)) * (3*x - 5)
    f_prime = function.diff(x)
    f_prime2 = f_prime.diff(x)
    f = lambdify(x, function)
    f_prime = lambdify(x, f_prime)
    f_prime2 = lambdify(x, f_prime2)
    a = 0
    b = 3
    epsilon = 0.00000001

    print("Run Newton - Raphson methode on F(x):")

    check_segments_f = (find_sign_change(split_segment(a, b), f))

    for n in check_segments_f:
        newton_raphson(n[0], n[1], f, f_prime, f, epsilon)

    print("Run Newton - Raphson methode on F'(x):")
    check_segments_f_prime = (find_sign_change(split_segment(a, b), f_prime))

    for n in check_segments_f_prime:
        newton_raphson(n[0], n[1], f_prime, f_prime2, f, epsilon)

    if f(0) == 0:
        print("The root is:", str(0) + '00000' + d + h + m)

    print("_____________________________________")
    print("Run Secant methode on F(x):")

    check_segments_f = (find_sign_change(split_segment(a, b), f))

    for n in check_segments_f:
        secant_method(n[0], n[1], f, f, epsilon)

    print("Run Secant methode on F'(x):")

    check_segments_f_prime = (find_sign_change(split_segment(a, b), f_prime))

    for n in check_segments_f_prime:
        secant_method(n[0], n[1], f_prime, f, epsilon)

    if f(0) == 0:
        print("The root is:", str(0) + '00000' + d + h + m)

    a = 0.5
    b = 1
    n = 30

    print("_____________________________________")
    print("Run Simpson methode:")
    print("Approximate area by Simpson", str(simpson(f, split_segment2(a, b, n), hight(a, b, n))) + '00000' + d + h + m)

    print("_____________________________________")
    print("Romberg methode:")
    print("Approximate area by Romberg", str(rombergMethod(function, a, b, 4, epsilon)) + '00000' + d + h + m)


main_methode()

