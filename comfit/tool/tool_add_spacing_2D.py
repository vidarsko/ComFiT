
def tool_add_spacing_2D(X,Y,U,V,spacing):
        X = X[::spacing, ::spacing]
        Y = Y[::spacing, ::spacing]
        U = U[::spacing, ::spacing]
        V = V[::spacing, ::spacing]
        return X,Y,U,V