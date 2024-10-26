module SpinTests

# Write your package code here.

end

using LinearAlgebra, StaticArrays, KomaMRI
function build_phantom()
    phantom = zeros(Bool, 100,100,100)
    r = 50
    for i in 1:100
        for j in 1:100
            for k in 1:100
                if (i-50)^2 + (j-50)^2 + (k-50)^2 < r^2
                    phantom[i,j,k] = 1
                end
            end
        end
    end
    return phantom
end

phantom = build_phantom()

off_resonance = zeros(100,100,100)


L = 5
basis = SolidHarmonics(L)
Rs = vec([SVector{3, Float64}(Tuple(x)) for x in CartesianIndices(off_resonance)])
Z = basis(Rs)

sph_basis= reshape(Z, (size(off_resonance)..., 36))
sph_field = zeros(100,100,100,36)
coeffs = randn(36)
function mul_sph_coeffs(sph_basis, coeffs, sph_field)
    for i in 1:36
       @views sph_field[:,:,:,i] .= sph_basis[:,:,:,i] .* coeffs[i] 
    end
end

mul_sph_coeffs(sph_basis, coeffs, sph_field)

off_resonance = sum(sph_field, dims=4)
scale = 1000
Rs_array = Array
rs = [getindex.(Rs, i)./1000 for i in 1:3]
pind = vec(phantom)
spins = sum(pind)
obj = Phantom{Float64}(x=rs[1][pind],y=rs[2][pind],z=rs[3][pind], T1=[100e-3 for _ in 1:spins], T2=[100e-3], Δw=vec(off_resonance)[pind], ρ=ones(length(pind))[pind])

#hyperfine scanner parameters
scanner = Scanner(0.0633, 1.0e-6, 20.0e-3, 200, 2.0e-6, 1.0e-5, 1.0e-5, 1.0e-6, 2.0e-5, 1e-4, 1e-5)

sys = Scanner(0.0633, 1.0e-6, 20.0e-3, 200, 2.0e-6, 1.0e-5, 1.0e-5, 1.0e-6, 2.0e-5, 1e-4, 1e-5)
ampRF = 2e-6                        # 2 uT RF amplitude
durRF = π / 2 / (2π * γ * ampRF)    # required duration for a 90 deg RF pulse
durRefocus = π / (2π * γ * ampRF) 
exc = [RF(ampRF,durRF);;]
refocus = [RF(ampRF,durRefocus);;]
ref_grad = [Grad(0.0, durRefocus);;]
Inv_grad = [Grad(0.0, durRF);;]
Inversion = Sequence(Inv_grad, exc)
Delay(t) = Sequence([Grad(0.0, t);;])

readout = [Grad(10e-3, 1e-3);;]

"""
    k_gradient(k_vector)

TBW
"""
function k_gradient(k_vector, Scanner)
    kx, ky, kz = k_vector
    Gmax = Scanner.Gmax
    Smax = Scanner.Smax
end

    

    



Inversion + Delay(1e-3)








function fast_spin_echo_module(FOV::Real, Dims=(100,100,100), sys)
    Δt = sys.ADC_Δt
	Gmax = sys.Gmax
    


end



function sequence_example(FOV::Real, N::Integer)

    # Define initial paramters (TODO: consider when N is even)
    sys = Scanner(0.0633, 1.0e-6, 2.0e-3, 200, 2.0e-6, 1.0e-5, 1.0e-5, 1.0e-6, 2.0e-5, 1e-4, 1e-5)
	Δt = sys.ADC_Δt
	Gmax = sys.Gmax
	Nx = Ny = N #Square acquisition
	Δx = FOV/(Nx-1)
	Ta = Δt*(Nx-1) #4-8 us
    Δτ = Ta/(Ny-1)
	Ga = 1/(γ*Δt*FOV)
	ζ = Ga / sys.Smax
	Ga ≥ sys.Gmax ? error("Ga=$(Ga*1e3) mT/m exceeds Gmax=$(Gmax*1e3) mT/m, increase Δt to at least Δt_min="*string(round(1/(γ*Gmax*FOV),digits=2))*" us.") : 0
	ϵ1 = Δτ/(Δτ+ζ)

	# EPI base
	EPI = Sequence(vcat(
	    [mod(i,2)==0 ? Grad(Ga*(-1)^(i/2),Ta,ζ) : Grad(0.,Δτ,ζ) for i=0:2*Ny-2],  #Gx
	 	[mod(i,2)==1 ? ϵ1*Grad(Ga,Δτ,ζ) :         Grad(0.,Ta,ζ) for i=0:2*Ny-2])) #Gy
	EPI.ADC = [mod(i,2)==1 ? ADC(0,Δτ,ζ) : ADC(N,Ta,ζ) for i=0:2*Ny-2]

	# Pre-wind and wind gradients
	ϵ2 = Ta/(Ta+ζ)
    PHASE =   Sequence(reshape(1/2*[Grad(      -Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) # This needs to be calculated differently
	DEPHASE = Sequence(reshape(1/2*[Grad((-1)^N*Ga, Ta, ζ); ϵ2*Grad(-Ga, Ta, ζ)],:,1)) # for even N
	seq = PHASE + EPI + DEPHASE

	# Saving parameters
	seq.DEF = Dict("Nx"=>Nx,"Ny"=>Ny,"Nz"=>1,"Name"=>"epi")

    # Return the sequence
	return seq
end

seq = sequence_example(0.2, 100)
plot_seq(seq; range=[0 30])
plot_phantom_map(obj, :ρ)

