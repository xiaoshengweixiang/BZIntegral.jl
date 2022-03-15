module BZInt3D
using DoubleFloats
using StaticArrays
using Base.Threads
export Lin3DRuleΘ,Lin3DRuleδ,Quad3DRuleΘ,Quad3DRuleδ,
       Quad3DRuleΘ𝔇,Quad3DRuleΘδ,Quad3DRuleΘΘ,Quad3DRuleΘΘ𝔇,
       Quad3DRuleδδ,Quad3DRuleΘΘδ,Quad3DRule𝒲,Quad3DRule𝒲𝒲,Quad3DRule𝒲𝒲𝒲,Quad3DRule𝒲𝔇,Quad3DRule𝒲𝒲𝔇
       
       
include("SplitMesh.jl")
include("QuadTetra.jl")
include("TetraSupply.jl")

"""
Linear tetrahedron method for weight function W(k) = Θ(eF-E(k))
with Blöchl correction
"""
function Lin3DRuleΘ(Emesh,eF)
    Tetras = Mesh2Tetra(size(Emesh)...)
    ETetras = Emesh[Tetras]
    WTetras = zeros(size(ETetras)...)

    @views @threads for i in 1:size(Tetras,1)
        WTetras[i,:] = LinTetraΘ_Blöchl(SVector{4}(ETetras[i,:]),eF)
    end

    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(Tetras,1)
        Wmesh[Tetras[i,:]] += WTetras[i,:]
    end

    return Wmesh/size(Tetras,1)*6
end
"""
Linear tetrahedron method for weight function W(k) = δ(eF-E(k))
"""
function Lin3DRuleδ(Emesh,eF)
    Tetras = Mesh2Tetra(size(Emesh)...)
    ETetras = Emesh[Tetras]
    WTetras = zeros(size(ETetras)...)

    @views @threads for i in 1:size(Tetras,1)
        WTetras[i,:] = LinTetraδ(SVector{4}(ETetras[i,:]),eF)
    end

    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(Tetras,1)
        Wmesh[Tetras[i,:]] += WTetras[i,:]
    end

    return Wmesh/size(Tetras,1)*6
end
"""
Recursive tetrahedron method for weight function W(k) = Θ(eF-E(k))
"""
function Quad3DRuleΘ(Emesh,eF,iter=2)
    QTetras = Mesh2QTetra(size(Emesh)...)
    EQTetras = Emesh[QTetras] 
    WQTetras = zeros(size(EQTetras)...) 
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraΘ(SVector{10}(EQTetras[i,:]),eF,iter)
    end

    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end
"""
Recursive tetrahedron method for weight function W(k) = W(k) = 1/D(k) Θ(eF-E(k))
"""
function Quad3DRuleΘ𝔇(Emesh,eF,Dmesh,iter=2)
    QTetras = Mesh2QTetra(size(Emesh)...)
    EQTetras = Emesh[QTetras]
    DQTetras = Dmesh[QTetras]
    WQTetras = zeros(size(EQTetras)...)
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraΘ𝔇(SVector{10}(EQTetras[i,:]),eF,SVector{10}(DQTetras[i,:]),iter)
    end
    
    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end
"""
Recursive tetrahedron method for weight function W(k) = δ(eF-E(k))
"""
function Quad3DRuleδ(Emesh,Kxmesh,Kymesh,Kzmesh,eF,iter=2)
    QTetras = Mesh2QTetra(size(Emesh)...)
    EQTetras = Emesh[QTetras]
    KxTetras = Kxmesh[QTetras]
    KyTetras = Kymesh[QTetras]
    KzTetras = Kzmesh[QTetras]
    for i in 1:size(QTetras,1)
       effm[i,:] = tetrainsert(KXTetra[i,:],KYTetra[i,:],KZTetra[i,:])
    WQTetras = zeros(size(EQTetras)...)
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraδ(SVector{10}(EQTetras[i,:]),eF,iter)
    end
    for i in 1:size(QTetras,1)
        for j in 1:size(WQTetra,2)
             weight_mass = weight+WQTetras[i,j]
        end
        total_weight = total_weight + weight_mass
        mbar = mabar + effmi[i,:]*weight_mass
    end
    return mbar/total_weight
   """ Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6"""
end

function tetrainsert(kxmesh,kymesh,kzmesh,emesh)
    b[]={0.0 for i in 1:10}
    a[][]={{0.0for i in 1:10} for i in 1:10}
    n=10
    for i in 1:10
         a[i][1]=kxmesh[i]*kxmesh[i]
         a[i][2]=2*kxmesh[i]*kymesh[i]
         a[i][3]=2*kxmesh[i]*kzmesh[i]
         a[i][4]=kymesh[i]*kymesh[i]
         a[i][5]=2*kymesh[i]*kzmesh[i]
         a[i][6]=kzmesh[i]*kzmesh[i]
         a[i][7]=kxmesh[i]
         a[i][8]=kymesh[i]
         a[i][9]=kzmesh[i]
         a[i][10]=1
         b[i]=emesh[i]
     end
     for k in 1:10
         max=abs(a[k][k])
         m=k
         for i in k+1:n-1
            if max < abs(a[i][k]):
                max = abs(a[i][k])
                m = i
            end
         end
         if m != k:
            for j in k:n-1
                t = a[m][j]
                a[m][j] = a[k][j]
                a[k][j] = t
            end
            t = b[m]
            b[m] = b[k]
            b[k] = t
         end    
        for i in k + 1:n-1
            for j in k+1:n-1
                a[i][j] -= a[k][j] * a[i][k] / a[k][k]
            end
            b[i] -= b[k] * a[i][k] / a[k][k]
            a[i][k] = 0
        end
    end
    x[n - 1] = b[n - 1] / a[n - 1][n - 1]

    for i in n-2:-1:0
        sum = 0.0
        for j in i + 1:n-1
            sum += a[i][j] * x[j]
        end
        x[i] = (b[i] - sum) / a[i][i]
    end
    print(x)
    return x
end
"""
Recursive tetrahedron method for weight function W(k) = δ(D(k)) Θ(eF-E(k))
"""
function Quad3DRuleΘδ(Emesh,eF,Dmesh,iter=2)
    QTetras = Mesh2QTetra(size(Emesh)...)
    EQTetras = Emesh[QTetras]
    DQTetras = Dmesh[QTetras]
    WQTetras = zeros(size(EQTetras)...)
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraΘδ(SVector{10}(EQTetras[i,:]),eF,SVector{10}(DQTetras[i,:]),iter)
    end
    
    Wmesh = zeros(typeof(float(eF)),size(Emesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = Θ(x1(k))*Θ(x2(k))
"""
function Quad3DRuleΘΘ(X1mesh,X2mesh,iter=2)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras]
    X2QTetras = X2mesh[QTetras]
    WQTetras = zeros(size(X1QTetras)...)
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraΘΘ(SVector{10}(X1QTetras[i,:]),SVector{10}(X2QTetras[i,:]),iter)
    end

    Wmesh = zeros(typeof(float(X1mesh[1])),size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = 1/D(k) Θ(x1(k))*Θ(x2(k))
"""
function Quad3DRuleΘΘ𝔇(X1mesh,X2mesh,Dmesh,iter=2)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras]
    X2QTetras = X2mesh[QTetras]
    DQTetras = Dmesh[QTetras]
    WQTetras = zeros(size(X1QTetras)...)
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraΘΘ𝔇(SVector{10}(X1QTetras[i,:]),SVector{10}(X2QTetras[i,:]),SVector{10}(DQTetras[i,:]),iter)
    end

    Wmesh = zeros(typeof(float(X1mesh[1])),size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = δ(D(k)) Θ(x1(k))*Θ(x2(k))
"""
function Quad3DRuleΘΘδ(X1mesh,X2mesh,Dmesh,iter=2)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras]
    X2QTetras = X2mesh[QTetras]
    DQTetras = Dmesh[QTetras]
    WQTetras = zeros(size(X1QTetras)...)
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetraΘΘδ(SVector{10}(X1QTetras[i,:]),SVector{10}(X2QTetras[i,:]),SVector{10}(DQTetras[i,:]),iter)
    end

    Wmesh = zeros(typeof(float(X1mesh[1])),size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end
"""
Recursive tetrahedron method for weight function W(k) = δ(x1(k))*δ(x2(k))
"""
function Quad3DRuleδδ(X1mesh,X2mesh,iter=2)
    dx1 = maximum(abs.(X1mesh))*5.0e-4
    xF = zero(X1mesh[1])
    Wmesh0 = Quad3DRuleΘδ(X1mesh,xF,X2mesh,iter)
    Wmeshd = Quad3DRuleΘδ(X1mesh,dx1+xF,X2mesh,iter)
    Wmesh = (Wmeshd-Wmesh0)/dx1
    return Wmesh
end


"""
Recursive tetrahedron method for weight function W(k) = 𝒲(x1(k))
default integral type is wtype=Float64
"""
function Quad3DRule𝒲(𝒲,X1mesh,iter=2,wtype=Float64)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras] 
    WQTetras = zeros(size(X1QTetras)...) 
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetra𝒲(𝒲,SVector{10}(X1QTetras[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))
default integral type is wtype=Float64
"""
function Quad3DRule𝒲𝒲(𝒲1,𝒲2,X1mesh,X2mesh,iter=2,wtype=Float64)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras] 
    X2QTetras = X2mesh[QTetras] 
    WQTetras = zeros(size(X1QTetras)...) 
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetra𝒲𝒲(𝒲1,𝒲2,SVector{10}(X1QTetras[i,:]),SVector{10}(X2QTetras[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))*𝒲3(x3(k))
default integral type is wtype=Float64
"""
function Quad3DRule𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,X1mesh,X2mesh,X3mesh,iter=2,wtype=Float64)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras] 
    X2QTetras = X2mesh[QTetras] 
    X3QTetras = X3mesh[QTetras] 
    WQTetras = zeros(size(X1QTetras)...) 
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetra𝒲𝒲𝒲(𝒲1,𝒲2,𝒲3,SVector{10}(X1QTetras[i,:]),SVector{10}(X2QTetras[i,:]),SVector{10}(X3QTetras[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = 𝒲(x1(k))/D(k)
default integral type is wtype=Float64
"""
function Quad3DRule𝒲𝔇(𝒲,X1mesh,Dmesh,iter=2,wtype=Float64)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras] 
    DQTetras = Dmesh[QTetras] 
    WQTetras = zeros(size(X1QTetras)...) 
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetra𝒲𝔇(𝒲,SVector{10}(X1QTetras[i,:]),SVector{10}(DQTetras[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

"""
Recursive tetrahedron method for weight function W(k) = 𝒲1(x1(k))*𝒲2(x2(k))/D(k)
default integral type is wtype=Float64
"""
function Quad3DRule𝒲𝒲𝔇(𝒲1,𝒲2,X1mesh,X2mesh,Dmesh,iter=2,wtype=Float64)
    QTetras = Mesh2QTetra(size(X1mesh)...)
    X1QTetras = X1mesh[QTetras] 
    X2QTetras = X2mesh[QTetras] 
    DQTetras = Dmesh[QTetras] 
    WQTetras = zeros(size(X1QTetras)...) 
    @views @threads for i in 1:size(QTetras,1)
        WQTetras[i,:] = QuadTetra𝒲𝒲𝔇(𝒲1,𝒲2,SVector{10}(X1QTetras[i,:]),SVector{10}(X2QTetras[i,:]),SVector{10}(DQTetras[i,:]),iter,wtype)
    end

    Wmesh = zeros(size(X1mesh)...)
    @views for i in 1:size(QTetras,1)
        Wmesh[QTetras[i,:]] += WQTetras[i,:]
    end
    return Wmesh/size(QTetras,1)*6
end

function fϵk(ϵ_μ,β)
    x = ϵ_μ*β
    if x>20 
        res = 0.0
    elseif x<-20
        res = 1.0
    else
        res = 1.0/(1+exp(x))
    end
    return res
end

function dfϵk_dϵ(ϵ_μ,β)
    x = ϵ_μ*β
    if x>20 || x<-20 
        res = 0.0
    else
        e =exp(x) 
        res = -β*e/(1+e)^2
    end
    return res
end
 
# Quad3DRule𝑓(Emesh,μ,β,iter=2)=Quad3DRule𝒲(e->fϵk(e-μ,β),Emesh,iter)

# Quad3DRule𝑑𝑓(Emesh,μ,β,iter=2)=Quad3DRule𝒲(e->dfϵk_dϵ(e-μ,β),Emesh,iter)

# Quad3DRule𝑑𝑓𝑑𝑓(E1mesh,E2mesh,μ,β,iter=2)=Quad3DRule𝒲𝒲(e->dfϵk_dϵ(e-μ,β),e->dfϵk_dϵ(e-μ,β),E1mesh,E2mesh,iter)

# Quad3DRule𝑓𝔓(Emesh,Dmesh,μ,β,η,iter=2)=Quad3DRule𝒲𝒲(e->dfϵk_dϵ(e-μ,β),d->1/(d+1im*η),Emesh,Dmesh,iter,ComplexF64)

# Quad3DRule𝑓𝑓𝔓(E1mesh,E2mesh,Dmesh,μ,β,η,iter=2)=Quad3DRule𝒲𝒲𝒲(e->dfϵk_dϵ(e-μ,β),e->dfϵk_dϵ(e-μ,β),d->1/(d+1im*η),E1mesh,E2mesh,Dmesh,iter,ComplexF64)


end

