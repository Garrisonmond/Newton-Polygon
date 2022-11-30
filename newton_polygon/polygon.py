import matplotlib.pyplot as plt
import seaborn as sns
from wolframclient.evaluation import WolframLanguageSession


class Fraction:
    """!
    Basic fraction class.
    """

    def __init__(self, a: int, b=1):
        """!
        Basic fraction class initializer
        @param a int enumerator
        @param b int denominator
        """
        if type(a) != int or type(b) != int:
            raise ValueError("The given fraction is improper")
        if b == 0:
            raise ValueError("Denominator is equal to 0")
        if b < 0:
            self.a = -a
            self.b = -b
        else:
            self.a = a
            self.b = b

    def __str__(self):
        """!
        Basic convert to string method.
        @return res str
        """
        if self.a == 0:
            return str(0)
        if self.b == 1:
            return str(self.a)
        return f'{self.a}/{self.b}'

    def __add__(self, other):
        """!
        Basic add method
        @other Fraction
        @return res Fraction
        """
        return Fraction(self.a * other.b + self.b * other.a, self.b * other.b)

    def __sub__(self, other):
        """!
        Basic sub method
        @other Fraction
        @return res Fraction
        """
        return Fraction(self.a * other.b - self.b * other.a, self.b * other.b)

    def __mul__(self, other):
        """!
        Basic mul method
        @other Fraction
        @return res Fraction
        """
        return Fraction(self.a * other.a, self.b * other.b)

    def __float__(self):
        """!
        Basic conver to float method
        @return res float
        """
        return self.a / self.b


def alphabetic(string: str):
    for letter in string:
        if not letter.isalpha():
            return 0
    return 1


class Factor:
    """!
    Class for handling basic factors inside a monomial.
    """

    def __init__(self, name: str, power: Fraction, der_ord=0, der_by=None):
        """!
        Factor class initializer.
        @param name str
        @param power Fraction
        @param der_ord int
        @param der_by str
        """
        if der_ord < 0 or type(der_ord) != int:
            raise ValueError("Passed Derivative order is not valid")
        self.der_ord = der_ord
        if type(name) != str:
            raise ValueError("Passed name is not a string")
        self.name = name
        if der_by is None and der_ord != 0:
            raise ValueError("Missed derived by parameter")
        if der_by is not None:
            if type(der_by) != str:
                raise ValueError("Passed derived by is not a string")
            if der_by == name:
                raise ValueError("The function is derived by itself")
            self.der_by = der_by
        else:
            self.der_by = ''
        if type(power) not in [int, Fraction]:
            raise ValueError("Passed power is invalid")
        if type(power) != Fraction:
            self.power = Fraction(power)
        else:
            self.power = power

    def __str__(self):
        """!
        Basic convert to string method.
        @return string str
        """
        string = self.name
        if self.power.a == 0:
            return ' '
        if self.der_ord != 0:
            for i in range(self.der_ord):
                string += "'"
        if self.der_by != '':
            string += f"[{self.der_by}]"
        if self.power.a != 1 or self.power.b != 1:
            string += f"^({str(self.power)})"

        return string


class Monomial:
    """!
    Class for handling differential monomials.
    """

    def __init__(self, factor_list: list):
        """!
        Monomial class initializer.
        @param factor_list list of Factor objects
        """
        for i in range(len(factor_list)):
            if type(factor_list[i]) != Factor:
                raise ValueError("Passed factor list is not valid")
        self.factor_list = sorted(sorted(
            factor_list, key=lambda x: x.name, reverse=True), key=lambda x: x.der_ord, reverse=True)

    def __str__(self):
        """!
        Conver to string method.
        @return string str
        """
        string = ''
        for factor in self.factor_list:
            s = factor.__str__()
            if s != ' ':
                string += ' '
                string += s

        return string

    def get_point(self, func: str, arg: str):
        """!
        Get a point used to initialize Newton Polygon. It is import to initialize a function and an argument,
        which is used to take derivative.
        @param func str
        @param arg str
        return [x, y] [Fraction, Fraction]
        """
        if type(func) != str:
            raise ValueError("Passed function is not a string")
        if not alphabetic(func):
            raise ValueError("Passed function is not alphabetic")
        if type(arg) != str:
            raise ValueError("Passed argument is not a string")
        if not alphabetic(func):
            raise ValueError("Passed argument is not alphabetic")
        if arg == func:
            raise ValueError("Passed argument is equal to a function")

        x = Fraction(0)
        y = Fraction(0)
        for factor in self.factor_list:
            if factor.name == func:
                if factor.der_by == arg:
                    x = x - Fraction(factor.der_ord) * factor.power
                y = y + factor.power
            if factor.name == arg and factor.der_by == '':
                x = x + factor.power
        return (x, y)


class WolframExpression:
    """!
    Class for basic for handling basic Wolfram expressions.
    """

    def __init__(self, wolfram_expr: str):
        """!
        WolframExpression initializer.
        @param wolfram_expr str
        """
        self.__session = WolframLanguageSession()
        if type(wolfram_expr) != str:
            ValueError("The given wolfram expression is not a string.")
        self.wolfram_expr = str(self.__session.evaluate(wolfram_expr))
        if self.wolfram_expr == '':
            ValueError("The given wolfram expression is empty.")
        self.__session.terminate()

        self.funcs = []
        self.args = []

    @staticmethod
    def _plus_str(args):
        """!
        Method for representing wolfram Plus[] as a string.
        @params list of str objects
        @return res str
        """
        if not len(args):
            return ''
        res = '('
        for i in range(0, len(args) - 1):
            res += args[i] + ' + '
        res += args[-1]
        res += ')'
        return res

    @staticmethod
    def _times_str(args: list):
        """!
        Method for representing wolfram Times[] as a string.
        @params list of str objects
        @return res str
        """
        res = ''
        if not len(args):
            return ''
        for i in range(0, len(args) - 1):
            res += args[i] + ' '
        res += args[-1]
        return res

    @staticmethod
    def _derivative_str(args: list):
        """!
        Method for representing wolfram Derivative[] as a string.
        @params list of str objects
        @return res str
        """
        if len(args) < 3:
            return ''
        return args[1] + "'" * int(args[0]) + f'[{args[2]}]'

    @staticmethod
    def _power_str(args: list):
        """!
        Method for representing wolfram Power[] as a string.
        @params list of str objects
        @return res str
        """
        if len(args) < 2:
            return ''
        return f'{args[0]}^({args[1]})'

    @staticmethod
    def _rational(args: list):
        """!
        Method for representing wolfram Rational[] as a string.
        @params list of str objects
        @return res str
        """
        if len(args) < 2:
            return ''
        return Fraction(int(args[0]), int(args[1])).__str__()

    @staticmethod
    def _other(f, args):
        """!
        Method for representing other wolfram functions as a string.
        @params list of str objects
        @return res str
        """
        if not len(args):
            return ''
        return f'{f}[{args[0]}]'


class Polynomial(WolframExpression):
    """!
       Basic class for handling differential polynomials initialized by wolfram expression.
    """
    def __init__(self, wolfram_expr: str):
        """!
        Polynomial initializer
        @param wolfram_expr str
        """
        super().__init__(wolfram_expr)
        self.monomial_list = []
        self.__factor_list = []
        self._py_expr, self.args = self._get(self.wolfram_expr)

    @staticmethod
    def _parse_factor(factor: str):
        """!
        Protected method used to parse factor expression.
        @param factor str
        @return name str
        @return power Fraction
        @return der_ord int
        @return der_by str
        """
        if type(factor) != str:
            ValueError("Passed factor is not a string.")
        name = ''
        power_a = '1'
        power_b = '1'
        der_by = ''
        der_ord = 0
        flag_name = 1
        flag_der_by = 0
        flag_power = 0
        flag_b = 0
        for s in factor:
            if s in ['[', ']', '^', "'", '(', ')', '/']:
                flag_name = 0
                if s == "'":
                    der_ord += 1
                elif s == '[':
                    flag_der_by = 1
                    continue
                elif s == ']':
                    flag_der_by = 0
                    continue
                elif s == '^':
                    power_a = ''
                    flag_power = 1
                    continue
                elif s == '/' and flag_power:
                    power_b = ''
                    flab_b = 1
                    continue
                else:
                    continue
            if flag_power:
                if not flag_b:
                    power_a += s
                else:
                    power_b += s
            if flag_der_by:
                der_by += s
            if flag_name:
                name += s
        power = Fraction(int(power_a), int(power_b))
        return name, power, der_ord, der_by

    def _factorize(self, monomial: str):
        """!
        Method used to factorize given monomial. Appends
        @param monomial str
        """

        for factor in monomial.split(' '):
            if type(monomial) != str:
                ValueError("Passed monomial is not a string.")
            name, power, der_ord, der_by = self._parse_factor(factor)
            self.__factor_list.append(Factor(name, power, der_ord, der_by))

    def _monomialize(self, monomials: list):
        """!
        Method used to convert monomial strings to Monomial objects.
        @param monomials list of str objects
        """
        for monomial in monomials:
            self._factorize(monomial)
            self.monomial_list.append(Monomial(self.__factor_list))
            self.__factor_list = []

    def _use_function(self, args: list, f):
        """!
        Function that calls initialized WolframExpression functions based on given name, and converts the output to str.
        @param args list of str
        @param f str
        @return str
        """
        for arg in args:
            if type(arg) != str:
                ValueError("Passed argument list is not valid.")
        if type(f) != str:
            ValueError("Passed function name is not a string.")
        if f == 'Plus':
            self._monomialize(args)
            return self._plus_str(args)
        if f == 'Derivative':
            return self._derivative_str(args)
        if f == 'Times':
            return self._times_str(args)
        if f == 'Power':
            return self._power_str(args)
        if f == 'Rational':
            return self._rational(args)
        return self._other(f, args)

    @staticmethod
    def _parse(wolfram_expr):
        """!
        Method used to parse initial wolfram expression.
        @param wolfram_expr str
        @return args list of str
        @return f str
        """
        left = 0
        args = []
        name = ''
        counter = 0
        f = ''
        for s in wolfram_expr:
            if name in ["Global`", "[Global`"]:
                name = ''

            if s == '[':
                left += 1
                if left == 1:
                    if f == 'Derivative':
                        left += 1
                    continue
            if s == ']':
                left -= 1
                if f == 'Derivative' and left == 1:
                    args.append(name)
                    name = ''
                    continue
                if left == 0:
                    args.append(name)
                    break

            if not left:
                f += s
                continue

            if left == 1 and s == ' ':
                continue

            if s == ',' and left == 1:
                args.append(name)
                counter += 1
                name = ''
                continue

            name += s
        return args, f

    def _get(self, wolfram_expr):
        """!
        Method used to recursively parse initial wolfram expression.
        @param wolfram_expr str
        @return py_expr str
        @return args list of str
        """
        args, f = self._parse(wolfram_expr)
        py_expr = [0] * len(args)
        if not len(args):
            return f, f
        else:
            self.funcs.append(f)
            for i in range(len(args)):
                py_expr[i], args[i] = self._get(args[i])
        py_expr = self._use_function(py_expr, f)

        return py_expr, args

    @staticmethod
    def _alter_py_expr(string):
        """!
        Method that alters py_expr of polynomial and makes it more user friendly.
        @param string str
        @return better str
        """
        if type(string) != str:
            ValueError("Passed argument list is not valid.")
        better = ''
        it = iter(range(len(string) - 2))
        for i in it:
            if string[i: i + 3] == '+ -':
                better += '- '
                next(it)
                next(it)
                continue
            better += string[i]
        better += string[-2:]
        if better[0] == '(' and better[-1] == ')':
            better = better[1:-1]
        if better[0] == '-':
            better = '- ' + better[1:]
        return better

    def __str__(self):
        """!
        Converts polynomial to string.
        @return res str
        """
        return self._alter_py_expr(self._py_expr)

    def replacement(self, wolfram_expr):
        """!
        Allows user to input replacement using wolfram language function definition of u[t]:= y[t] + c/t^-1
        @param wolfram_expr str
        @return res Polynomial
        """
        session = WolframLanguageSession()
        session.evaluate(wolfram_expr)
        new = session.evaluate(f'Expand[{self.wolfram_expr}]')
        session.terminate()
        return Polynomial(new)


class NewtonPolygon(Polynomial):
    """!
    Class for handling basic Newton Polygons for differential polynomials.
    """
    def __init__(self, wolfram_expr: str, func: str, arg: str):
        """!
        Newton Polygon initializer.
        @param func str
        @param arg str
        @param wolfram_expr str
        @return res NewtonPolygon
        """
        super().__init__(wolfram_expr)
        self.points = []
        if type(func) != str:
            ValueError("Passed function name is not a string.")
        if type(arg) != str:
            ValueError("Passed argument name is not a string.")
        self.func = func
        self.arg = arg
        self.edges = []
        self._points = {}
        for monomial in self.monomial_list:
            point = monomial.get_point(func, arg)
            for coord in point:
                if type(coord) != Fraction:
                    ValueError("Passed points list is not valid.")

            if point in self._points.keys():
                self._points[point].append(monomial)
            else:
                self._points[point] = [monomial]

    def replacement(self, wolfram_expr):
        """!
        Allows user to input replacement using wolfram language function definition of u[t]:= y[t] + c/t^-1
        @param wolfram_expr str
        @return res Polynomial
        """
        session = WolframLanguageSession()
        session.evaluate(wolfram_expr)
        new = session.evaluate(f'Expand[{self.wolfram_expr}]')
        session.terminate()
        return NewtonPolygon(new, self.func, self.arg)

    def draw(self, name=''):
        """!
        Method for drawing Newton polygons.
        @param name str"""
        if type(name) != Fraction:
            ValueError("Passed name is not a string.")
        sns.set(style="darkgrid")
        l = list(self._points.keys())
        print(l)
        for i in range(len(l)):
            plt.plot(float(l[i][0]), float(l[i][1]), marker='o', label=f'Q{i}')
            plt.text(float(l[i][0]), float(l[i][1]),
                     f'({float(l[i][0])},{float(l[i][1])})')
        for i in range(len(self.edges)):
            point1 = self.points[self.edges[i][0]]
            point2 = self.points[self.edges[i][1]]
            x_values = [float(point1[0]), float(point2[0])]
            y_values = [float(point1[1]), float(point2[1])]
            plt.plot(x_values, y_values, label=f'Г1_{i}')

        plt.legend()
        plt.title(name)
        plt.show()

    def add_edge(self, point_index1: int, point_index2: int):
        """!
        Add edge to the Newton polygon.
        @param point_index1 int
        @param point_index2 int
        """
        if type(point_index1) != int or type(point_index2) != int:
            raise NameError("Passed indexes are not int")
        if point_index1 >= len(self.points) or point_index1 < 0:
            raise NameError("Passed index 1 is out of range")
        if point_index2 >= len(self.points) or point_index2 < 0:
            raise NameError("Passed index 2 is out of range")
        self.edges.append([point_index1, point_index2])

    def remove_edge(self, edge_index):
        """!
        Remove edge from the Newton polygon.
        @param edge_index int
        """
        if type(edge_index) != int:
            raise NameError("Passed index is not int")
        if edge_index >= len(self.edges) or edge_index < 0:
            raise NameError("Passed index is out of range")
        del self.edges[edge_index]

