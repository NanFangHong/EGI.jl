module EGI

using LinearAlgebra, Polyhedra, QHull, Makie
export ExtendedGaussImage, Surface, Volume, Polyhedron, Normal, randPolyhedron, RenderPolyhedron, Little, SavePolyhedronArray

struct ExtendedGaussImage
    Normal::Array{Float64,2} where T <: Union{Int,Float64}
    Area::Union{Array{T,2},Array{T,1}} where T <: Union{Int,Float64}
    function ExtendedGaussImage((Normal, Area))
        if size(Normal)[2] != length(Area)
            error("Number of surface normal does not match number of surface area. ")
        end
        if sum([norm(Normal[:,i]) for i in 1:length(Area)] - ones(length(Area))) > 1e-10
            error("Normal not unit vector")
        end
        if sum(Normal*Area) > 1e-10
            error("It is not a valid convex polyhedron (Does not satisfy sum of surfaces normal times area equals 0")
        end
        new(Normal,Area)
    end
end

function randEGI(faces::Int64)
    for i = 1:(10*faces)
        global Normed = zeros(3,1)
        for j = 1:faces
            N = rand(3)*2 .-1
            N = N ./ norm(N)
            Normed = [Normed N]
        end
        global Normed = Normed[:,2:end]

        v2 = rand(faces - 3, 1)
        v1 = inv(Normed[:, 1:3])*(-Normed[:, 4:end]*v2)
        global v = [v1;v2]
        if (sum(v .> zeros(faces,1)) == faces)
            break
        end
    end
    return Normed,v
end

function Surface(p::QHull.Polyhedron)
    S = Float64[]
    for pi in eachindex(halfspaces(p))
        global n = get(p, pi).a
        n = n/norm(n)
        get(p, pi).β
        global pts = incidentpoints(p, pi)
        global b = pts[2] - pts[1]
        b = b/norm(b)
        global a = cross(n, b)
        norm(a)
        append!(S, Float64(Polyhedra.volume(polyhedron(vrep([[a b]'*pt for pt in pts]), QHull.Library()))))
    end
    return S
end

function Volume(H::Array{Float64,1}, S::Array{Float64,1})
    return 1/3*(H'*S)
end

function Volume(p::QHull.Polyhedron)
    faces = length(hrep(p).b)
    H = MixedMatHRep(hrep(p)).b ./[norm(MixedMatHRep(hrep(p)).A[i,:]) for i in 1:faces]
    S = Float64[]
    for pi in eachindex(halfspaces(p))
        global n = get(p, pi).a
        n = n/norm(n)
        get(p, pi).β
        global pts = incidentpoints(p, pi)
        global b = pts[2] - pts[1]
        b = b/norm(b)
        global a = cross(n, b)
        norm(a)
        append!(S, Float64(Polyhedra.volume(polyhedron(vrep([[a b]'*pt for pt in pts]), QHull.Library()))))
    end
    return 1/3*(H'*S)
end

function Polyhedron(H,N)
    h = hrep(N', [norm(N[:,i]) for i in 1:length(H)] .* H)
    p = polyhedron(h, QHull.Library())
    return p
end

function Normal(p::QHull.Polyhedron)
    faces = length(hrep(p).b)
    return [MixedMatHRep(hrep(p)).A[i,:]./norm(MixedMatHRep(hrep(p)).A[i,:]) for i in 1:faces]
end

function randPolyhedron(faces)
    egi = randEGI(faces)
    img = ExtendedGaussImage(egi)
    for i in 1:(10*faces)
        H = rand(faces)
        global p = Polyhedron(H, img.Normal)
        if length(Surface(p)) == faces
            break
        end
    end
    return p, ExtendedGaussImage((img.Normal, Surface(p)[:,:]))
end

function RenderPolyhedron(p::QHull.Polyhedron)
    m = Polyhedra.Mesh(p)
    scene = mesh(m, color=:white)
end

Base.length(q::QHull.Polyhedron) = 1
Base.iterate(q::QHull.Polyhedron) = (q, nothing)
function Little(k, ϵ, img::ExtendedGaussImage)
    faces = length(img.Area)
    H = ones(faces)
    p1 = img.Area / sum(img.Area)
    KLArray = []
    KL = 0.05
    qArray = QHull.Polyhedron[]
    q = Polyhedron(H, img.Normal)
    while KL >= ϵ
        H = H / (Volume(q))^(1/3)
        q = Polyhedron(H, img.Normal)
        ∇ = 1/3*Surface(q)
        Δ = ∇ - (img.Area * norm(∇))/ norm(img.Area)
        H = vec(H + k .* Δ)
        H = H / (Volume(q))^(1/3)
        q = Polyhedron(H, img.Normal)
        p2 = Surface(q) / sum(Surface(q))
        KL = sum([p2[i]*log(p2[i]/p1[i]) for i in 1:faces])
        @show p2 ./ p1
        @show KL
        append!(KLArray, KL)
        append!(qArray, q)
    end
    return qArray, KLArray
end

function SavePolyhedronArray(qArray::Array{QHull.Polyhedron,1})
    for i in 1:length(qArray)
        scene = Scene()
        m = Polyhedra.Mesh(qArray[i])
        mesh!(scene, m, color=:white)
        Makie.save(string(i)*".png", scene)
        sleep(1/24)
    end
end

end # module
