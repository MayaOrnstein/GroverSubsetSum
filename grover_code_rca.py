#! /bin/env python
# Code for grover subset sum solver using a ripple carry adder (impmplented entirely with H, x, cx, ccx gates)

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister

import numpy as np


# Implemented based on https://arxiv.org/pdf/quant-ph/0410184
def unmajority_circuit(quantum_circuit, cpa, bpa, c):
    r'''
    $c_i$ The ith carry bit
    $b_i$ the ith bit of b
    $a_i$ the ith bit of a

    $c_i \oplus a_i \to c_i$
    $b_i \oplus a_i \oplus s_i$
    $c_{i+1} \to a_i$
    '''
    quantum_circuit.ccx(cpa, bpa, c)
    quantum_circuit.cx(c, cpa)
    quantum_circuit.cx(cpa, bpa)


# Majority
def maj_circuit():
    r'''
    .. math::
        $c_i$ The ith carry bit
        $b_i$ the ith bit of b
        $a_i$ the ith bit of a

        $c_i \to c_i \oplus a_i$
        $b_i \to b_i \oplus a_i$
        $a_i \to c_{i+1}$
    '''
    c_reg = QuantumRegister(size=1, name=r"c_i")
    b_reg = QuantumRegister(size=1, name=r"b_i")
    a_reg = QuantumRegister(size=1, name=r"a_i")
    maj = QuantumCircuit(c_reg, b_reg, a_reg, name=r"MAJ")
    maj.cx(a_reg[0], b_reg[0])
    maj.cx(a_reg[0], c_reg[0])
    maj.ccx(c_reg[0], b_reg[0], a_reg[0])
    return maj


# Unmajority and add
def uma_circuit(circuit=False):
    r'''
    $c_i$ The ith carry bit
    $b_i$ the ith bit of b
    $a_i$ the ith bit of a

    $c_i \oplus a_i \to c_i$
    $b_i \oplus a_i \oplus s_i$
    $c_{i+1} \to a_i$
    '''
    cpa_reg = QuantumRegister(size=1, name=r"c_i \oplus a")
    cpb_reg = QuantumRegister(size=1, name=r"c_i \oplus b")
    c_reg = QuantumRegister(size=1, name=r"c_\{i+1\}")

    uma = QuantumCircuit(cpa_reg, cpb_reg, c_reg, name="UMA")
    uma.ccx(cpa_reg[0], cpb_reg[0], c_reg[0])
    uma.cx(c_reg[0], cpa_reg[0])
    uma.cx(cpa_reg[0], cpb_reg[0])
    return uma


def ripple_carry_adder(n_qubits:int, overflow:bool=False):
    """
    A Quantum circuit that adds a to b (both n_qubits in length)

    |ab> -> |as> where s is the sum of a and b

    if overflow is set to True, we instead get
    |abz> -> |asz> where z is the overflow bit

    input wires:
    -
    - a -
    -

    -
    - b -
    -

    - c -
    """

    ripple_carry_add = QuantumCircuit(name="RCA")
    a_reg = QuantumRegister(size=n_qubits, name="a")
    b_reg = QuantumRegister(size=n_qubits, name="b")
    c_reg = AncillaRegister(size=1, name="carry_in")
    if overflow:
        z_reg = AncillaRegister(size=1, name="carry_out")
        ripple_carry_add.add_register(a_reg, b_reg, c_reg, z_reg)
    else:
        ripple_carry_add.add_register(a_reg, b_reg, c_reg)

    maj = maj_circuit().to_gate()
    uma = uma_circuit().to_gate()
    ripple_carry_add.compose(maj, qubits=[c_reg, b_reg[0], a_reg[0]], inplace=True)
    for i in range(1, n_qubits):
        ripple_carry_add.compose(maj, qubits=[a_reg[i-1], b_reg[i], a_reg[i]], inplace=True)
    if overflow:
        ripple_carry_add.cx(a_reg[n_qubits-1], z_reg[0])
    for i in range(n_qubits - 1, 0, -1):
        ripple_carry_add.compose(uma, qubits=[a_reg[i-1], b_reg[i], a_reg[i]], inplace=True)
    ripple_carry_add.compose(uma, qubits=[c_reg, b_reg[0], a_reg[0]], inplace=True)
    return ripple_carry_add


def mcz_circuit(n_qubits):
    r"""
    A multi-controlled z gate implemented using only Hadamard, cnot, not, and Tofelli gates

    input wires:
    -
    - ancilla register -
    -

    -
    - input -
    -

    - c (should be |0>)

    n_quibtis is the number of quibits in the mcz gate
    """
    # If we do not have enough qubits for a mcz circuit, do nothing
    if n_qubits <= 1:
        return QuantumCircuit(n_qubits*2 + 1, name="ID")

    mcz = QuantumCircuit(name="MCZ")
    a_reg = QuantumRegister(size=n_qubits-1, name="a")
    b_reg = QuantumRegister(size=n_qubits, name="b")
    c_reg = AncillaRegister(size=1, name="carry_in")

    mcz.add_register(a_reg, b_reg, c_reg)
    maj = maj_circuit().to_gate()
    maj_dag = maj.inverse()

    # Carry in a 1
    mcz.x(c_reg[0])
    mcz.compose(maj, qubits=[c_reg, b_reg[0], a_reg[0]], inplace=True)
    for i in range(1, n_qubits-1):
        mcz.compose(maj, qubits=[a_reg[i-1], b_reg[i], a_reg[i]], inplace=True)
    # If we overflowed most significant qubit we negate the final qubit
    #mcz.barrier()
    mcz.h(b_reg[n_qubits-1])
    mcz.cx(a_reg[n_qubits-2], b_reg[n_qubits-1])
    mcz.h(b_reg[n_qubits-1])
    #mcz.barrier()
    # Uncompute
    for i in range(n_qubits - 2, 0, -1):
        mcz.compose(maj_dag, qubits=[a_reg[i-1], b_reg[i], a_reg[i]], inplace=True)
    mcz.compose(maj_dag, qubits=[c_reg, b_reg[0], a_reg[0]], inplace=True)
    mcz.x(c_reg[0])
    return mcz


def controlled_add_circuit(n_qubits, val):
    """
    A controlled adder gate.

    input wires:
    - Control -

    -
    - ancilla registers -
    -

    -
    - running sum registers -
    -

    - Carry bit register (for the adder, presumably |0>) -


    n_qubits is the number of quibits for the running sum
    val is the value we are adding to the running sum
    """
    x_reg = QuantumRegister(size=1, name='x')
    sum_regs = QuantumRegister(size=n_qubits, name='s')
    anc_regs = AncillaRegister(size=n_qubits, name='a')
    c_reg = AncillaRegister(size=1, name='c')
    controlled_add = QuantumCircuit(
        x_reg,
        anc_regs,
        sum_regs,
        c_reg,
        name=f"Controlled Add {val}")

    # Twos-Complement representation of the value
    twos_complement = val if val >= 0 else 2**n_qubits + val
    bit_string = f'{int(twos_complement):>0{n_qubits}b}'

    rca = ripple_carry_adder(n_qubits, overflow=False).to_gate()

    # Set the value if x_i is 1
    for i, bit in enumerate(bit_string[::-1]):
        if bit == "1":
            controlled_add.cx(x_reg[0], anc_regs[i])

    # Add the input to the sum
    controlled_add.compose(rca, qubits=controlled_add.qubits[1:], inplace=True)

    # Reset the input
    for i, bit in enumerate(bit_string[::-1]):
        if bit == "1":
            controlled_add.cx(x_reg[0], anc_regs[i])

    return controlled_add


def weighted_sum_circuit(n_qubits, vals):
    """
    Computes the "weighted sum" (a sum of values (each adder is controlled by a control bit (i.e. the weight))
    input wires:
    -
    - Control registers -
    -

    -
    - ancilla registers -
    -

    -
    - running sum registers -
    -

    - c (|0>) -

    n_qubits is the number of quibits for the running sum
    vals is the list of values we are adding to the running sum (controled by corresponding control bits
    """
    n_x_vals = len(vals)
    x_regs = QuantumRegister(size=n_x_vals, name='x')
    s_regs = QuantumRegister(size=n_qubits, name='s')
    anc_regs = AncillaRegister(size=n_qubits, name='a')
    c_reg = AncillaRegister(size=1, name='c')

    weighted_sum = QuantumCircuit(
        x_regs,
        anc_regs,
        s_regs,
        c_reg,
        name=f"WS {vals}")

    for i, val in enumerate(vals):
        controlled_adder = controlled_add_circuit(n_qubits, val).to_gate()
        weighted_sum.compose(controlled_adder, [x_regs[i]] + weighted_sum.qubits[n_x_vals:], inplace=True)
    return weighted_sum


def subsetsum_oracle_circuit(n_qubits, vals, target):
    """
    Computes (weighted sum, reflection about the target vale, uncompute wighted sum)

    input wires:
    -
    - Control registers -
    -

    -
    - ancilla registers -
    -

    -
    - running sum registers (also ancilla) -
    -

    - c (|0>) -

    n_qubits is the number of quibits for the running sum
    vals is the list of values we are adding to the running sum (controled by corresponding control bits
    target is the target for the sum
    """

    n_x_vals = len(vals)
    x_regs = QuantumRegister(size=n_x_vals, name='x')
    s_regs = AncillaRegister(size=n_qubits, name='s')
    anc_regs = AncillaRegister(size=n_qubits, name='a')
    c_reg = AncillaRegister(size=1, name='c')

    subsetsum_oracle = QuantumCircuit(
        x_regs,
        anc_regs,
        s_regs,
        c_reg,
        name="SubsetSum Oracle")

    twos_complement = target if target >= 0 else 2**n_qubits + target
    bit_string = f'{int(twos_complement):>0{n_qubits}b}'

    weighted_sum = weighted_sum_circuit(n_qubits, vals).to_gate()
    weighted_sum_dag = weighted_sum.inverse()
    mcz = mcz_circuit(n_qubits).to_gate()

    subsetsum_oracle.compose(weighted_sum, inplace=True)
    # Set the value if x_i is 1
    for i, bit in enumerate(bit_string[::-1]):
        if bit == "0":
            subsetsum_oracle.x(s_regs[i])
    subsetsum_oracle.compose(mcz, qubits=(anc_regs[0:n_qubits-1] + s_regs[:] + c_reg[:]), inplace=True)
    for i, bit in enumerate(bit_string[::-1]):
        if bit == "0":
            subsetsum_oracle.x(s_regs[i])
    subsetsum_oracle.compose(weighted_sum_dag, inplace=True)

    return subsetsum_oracle


def prep_circuit(n_vals:int):
    """
    Prepares the control bits into |++...+>
    i.e.

    -- H --
    -- H --
    .
    .
    .
    -- H --

    n_vals is the number of control bits
    """
    quantum_register = QuantumRegister(size=n_vals, name='x')
    prep = QuantumCircuit(quantum_register, name="State Prep")
    prep.h(quantum_register)
    return prep


def grover_diffuser_circuit(n_vals):
    """
    Diffuser circuit for Grovers algorithm

    input wires:
    -
    - Control registers -
    -

    -
    - ancilla registers -
    -

    - c (|0>) -

    n_vals is the number of control bits
    """

    x_regs = QuantumRegister(size=n_vals, name='x')
    anc_regs = AncillaRegister(size=n_vals-1, name='a')
    c_reg = AncillaRegister(size=1, name='c')
    diffuser = QuantumCircuit(x_regs, anc_regs, c_reg, name="Diffuser")

    mcz = mcz_circuit(n_vals).to_gate()

    # Hadamard then Not gates
    diffuser.h(x_regs)
    diffuser.x(x_regs)
    # Multi-controlled z gate
    diffuser.compose(mcz, qubits=(anc_regs[:] + x_regs[:] + c_reg[:]), inplace=True)
    # Not gates then Haddamard gates
    diffuser.x(x_regs)
    diffuser.h(x_regs)

    return diffuser


def grover_subsetsum_circuit(vals, target, n_iter=1):
    """
    Full Circuit for the subset sum problem using grovers algorithm
    """
    n_qubits = int(
        np.ceil(
            max(
                np.log2(sum(np.abs(vals))),
                np.log2(np.abs(target))
            )
        ) + 1
    )
    n_vals = len(vals)
    n_ancilla = max(2*n_qubits, n_vals) + 1

    x_regs = QuantumRegister(size=n_vals, name='x')
    anc_regs = AncillaRegister(size=n_ancilla, name='a')

    prep = prep_circuit(n_vals).to_gate()
    oracle = subsetsum_oracle_circuit(n_qubits, vals, target).to_gate()
    diffuser = grover_diffuser_circuit(n_vals).to_gate()

    grover_subsetsum = QuantumCircuit(
        x_regs,
        anc_regs,
        name="Grover SubsetSum")

    grover_subsetsum.compose(prep, x_regs[:], inplace=True)
    for _ in range(n_iter):
        grover_subsetsum.compose(oracle, qubits=(x_regs[:] + anc_regs[:2*n_qubits+1]), inplace=True)
        grover_subsetsum.compose(diffuser, qubits=(x_regs[:] + anc_regs[:n_vals]), inplace=True)
    return grover_subsetsum
