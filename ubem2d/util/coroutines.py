__all__ = ['arithmetic_stepper', 'function_stepper']

def arithmetic_stepper(dx, x=0., stop=None):
    '''
    A coroutine which returns the result of incrementing a variable by a
    fixed amount on each call.  The value of the increment may be changed
    on the fly.
    '''
    step = 0
    while(True):
        if (stop is not None and stop(step, x)):
            raise StopIteration()
        dx_new = yield x
        if (dx_new is not None):
            dx = dx_new
        x += dx
        step += 1

def function_stepper(domain_generator, f):
    '''
    A generator which yields the graph of a given function over a domain.
    '''
    while (True):
        try:
            x = next(domain_generator)
            yield x, f(x)
        except:
            # https://stackoverflow.com/questions/51700960
            # https://www.python.org/dev/peps/pep-0479/
            return

if __name__ == '__main__':
    time = arithmetic_stepper(.1, stop = lambda i,x: i >= 10)
    f = lambda t: t**2
    for i, (t, ft) in enumerate(function_stepper(time, f)):
        print(('{:5d} ' + 2*'{:>10.4f}').format(i, t, ft))
