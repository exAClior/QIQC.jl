using MLStyle

@data Register begin
    # TODO needs to redesign
    BasicRegister(state::Int)
    CompositeRegister(sub_registers::AbstractVector{<:Register})
end

function size(r::Register)
    @match r begin
        BasicRegister(state) => 2^state
        CompositeRegister(sub_registers) => prod(size.(sub_registers))
    end
end

cr = CompositeRegister([BasicRegister(1), BasicRegister(2)])

size(cr)
