# from typing import NamedTuple
#
# class Abc(NamedTuple):
#     A: float
#     B: float
#     C: float
#
# Data = Abc(A = 1.0, B = 3.0, C = 5.0)
# A, B, C = Data
# A = 3.6
# print(A, ' ', B, ' ', C)
#
#
#
# from dataclasses import dataclass
#
# @dataclass
# class Abc():
#     A: float
#     B: float
#     C: float
#
# Data = Abc(A = 1.0, B = 3.0, C = 5.0)
# # A, B, C = Data  ----- так сделать не получится, т.к. возвращается не кортеж
# A = Data.A
# B = Data.B
# C = Data.C
# A = 3.6
# print(A, ' ', B, ' ', C)
# print(Data.A)



# class Person:
#     age: str
#     name: str
#
#
# Person_1 = Person()
# Person_1.age = '34'
# Person_1.name = 'John'
#
# Person_2 = Person()
# Person_2.age = '28'
# Person_2.name = 'Bob'
#
# print(Person_2.age)
#
#
class Person:
    def __init__(self, age, name):
        self.age = age
        self.name = name

    def Print(self):
        print("Name:", self.name, ";  Age:", self.age)

    def Print2(self):
        self.Print()
        print("Dreams: just wants be happy")


Person_1 = Person(age = '34', name = 'John')
# Person_1 = Person()
# Person_1.age = '34'
# Person_1.name = 'John'   ------   так сделать не получится из-за __init__

Person_2 = Person(age = '28', name = 'Bob')


Person_1.Print2()





















