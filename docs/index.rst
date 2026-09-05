=========================================
varietas — algebraic kinematics for ROS 2
=========================================

**Exact inverse kinematics, workspace implicitization, and singularity variety
decomposition via Gröbner bases and elimination theory.**

varietas treats robot kinematics as what it is — *a system of polynomial
equations over a finite-dimensional quotient ring* — and solves it with the
tools of computational algebraic geometry rather than with numerical iteration
or structural pattern matching.

Given a URDF, the library recovers the chain exactly over :math:`\Q`,
rationalizes the forward kinematics map, and computes a Gröbner basis of the
resulting ideal. Solutions are complete: the shape of the ideal certifies that
no branch has been missed, and the same ideal yields the implicit equations of
the reachable workspace and the singular locus as a variety with a dimension
and an image.

Where IKFast pattern-matches known kinematic structures and fails silently
outside them, varietas either produces a solution set with a proof of
completeness or reports precisely which hypothesis of the Extension Theorem was
violated.

.. figure:: figures/pipeline.svg
   :width: 100%
   :alt: The varietas pipeline, from a robot description to the three
         constructions the ideal supports

   **The pipeline.** A robot description is recovered exactly over the
   rationals, rationalised into an ideal, and reduced to a Gröbner basis; that
   one basis then supports all three constructions — the configurations, the
   workspace closure, and the singular locus. Everything up to the last step is
   exact, and the only numerical operation in the library is the
   eigendecomposition that reads the configurations off the action matrix.

.. rubric:: In one command

.. code-block:: sh

   # A URDF in, a self-contained C++ header out.
   ros2 run varietas_urdf urdf_codegen arm.urdf arm_ik.hpp

   # Or one pose, solved exactly, orientation included.
   ros2 run varietas_urdf urdf_solve arm.urdf --xyz 0.4 0.1 0.3 --rpy 0 1.57 0

.. rubric:: Where to start

:doc:`installation`
   Dependencies, ``colcon build``, and building without ROS at all.

:doc:`quickstart`
   From a URDF to a generated solver, and from a pose to a configuration.

:doc:`theory/index`
   What is actually being computed, and why the answers are complete.

:doc:`guide/index`
   The packages, the command line tools, and the shape of generated code.

:doc:`api/index`
   The headers, type by type.

:doc:`status`
   What works, what it costs, and what it refuses. Read this before
   planning around the library.

.. warning::

   varietas is at version 0.1.0 and its reach is deliberately narrow. The
   parametric path is limited to **two adjoined pose parameters**, and the
   fixed-pose path stops being interactive at **five joints**. A seven-axis
   manipulator is not a target this library can serve today; see
   :doc:`status` for exactly why and :doc:`roadmap` for what would change it.

.. toctree::
   :maxdepth: 2
   :caption: Getting started
   :hidden:

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Theory
   :hidden:

   theory/index

.. toctree::
   :maxdepth: 2
   :caption: User guide
   :hidden:

   guide/index

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   api/index

.. toctree::
   :maxdepth: 2
   :caption: Project
   :hidden:

   status
   roadmap
   changelog
   license

Indices
=======

* :ref:`genindex`
* :ref:`search`
