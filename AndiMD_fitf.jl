import LsqFit.curve_fit
using LsqFit
using ForwardDiff
"""
curve_fit where the model(x,p) function is transformed by with prescaled parameters p = scale_fwd(scale_bck(p)):
    model_scaled(x,p) = model(x, scale_bck(p))
    model_scaled(x,scale_fwd(p)) == model(x,p)
"""
function curve_fit(model, xdata::AbstractArray, ydata::AbstractArray, p0::AbstractArray, scale_fwd::Function, scale_bck::Function; inplace, kwargs...)
    
    model_scaled(x, p_fwd) = model(x, scale_bck(p_fwd))
    p0_fwd = scale_fwd(p0)
    res_fwd = curve_fit(model_scaled, xdata, ydata, p0_fwd; inplace=inplace, kwargs...)
    
    param_reconstructed = scale_bck(res_fwd.param)
    jac_fwd_inv = ForwardDiff.jacobian(scale_bck, res_fwd.param) |> inv
    resid_reconstructed = model.(xdata, (param_reconstructed,) ) .- ydata
    
    res_reconstructed = LsqFit.LsqFitResult(param_reconstructed, resid_reconstructed, res_fwd.jacobian*jac_fwd_inv, res_fwd.converged, res_fwd.wt)
end

