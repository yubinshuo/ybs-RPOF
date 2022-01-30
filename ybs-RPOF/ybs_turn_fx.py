import numpy as np

def func_transformer(func):
    prefered_function_format = '''
    def demo_func(x):
        x1, x2, x3 = x[:, 0], x[:, 1], x[:, 2]
        return x1 ** 2 + (x2 - 0.05) ** 2 + x3 ** 2
    '''
    is_vector = getattr(func, 'is_vector', False)
    if is_vector:
        return func
    else:
        if func.__code__.co_argcount == 1:
            def func_transformed(X):
                return np.array([func(x) for x in X])
            return func_transformed
        elif func.__code__.co_argcount > 1:
            def func_transformed(X):
                return np.array([func(*tuple(x)) for x in X])
            return func_transformed
    raise ValueError('''
    object function error,
    function should be like this:
    ''' + prefered_function_format)
