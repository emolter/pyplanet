class DataReturn:
    f = []
    b = []
    Tb = []
    header = {}

    def __repr__(self):
        s = ''
        for i, b in enumerate(self.b):
            bstr = 'b = {}:'.format(b)
            s += bstr + '\n  '
            f = ['{:6.1f}'.format(x) for x in self.f]
            s += ' '.join(f)
            s += '\n  '
            T = ['{:6.1f}'.format(x) for x in self.Tb[i]]
            s += ' '.join(T) + '\n'
        return s
