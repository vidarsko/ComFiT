

def tool_add_spacing_3D(X,Y,Z,U,V,W,spacing):
        X = X[::spacing, ::spacing, ::spacing]
        Y = Y[::spacing, ::spacing, ::spacing]
        Z = Z[::spacing, ::spacing, ::spacing]
        U = U[::spacing, ::spacing, ::spacing]
        V = V[::spacing, ::spacing, ::spacing]
        W = W[::spacing, ::spacing, ::spacing]
        return X,Y,Z,U,V,W  