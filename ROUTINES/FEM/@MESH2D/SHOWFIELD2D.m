function [] = SHOWFIELD2D(m, varargin)
%SHOWFIELD2D Depicts the field on the 2D mesh. All inputs are optional.
%
%  USAGE:
%  ------
%    MESH.SHOWFIELD2D(U);
%  INPUTS:
%  -------
%    U		: Nil (for just mesh)
%		  (MESH.Nn,1) for nodal field
%  		  (MESH.Ne,1) for element field
%  		  (MESH.Ne*MESH.Nq^2,1) for tractions at quadrature locations
%    e_on	: vector of indices of elements to plot.
%  OUTPUTS:
%  --------
%    Nil

  if isempty(varargin)
    U = 1;
  else
    U = varargin{1};
  end
  if length(varargin)<2
    e_on = 1:m.Ne;
  else
    e_on = varargin{2};
  end
  if size(U, 2)==1
      alpha = ones(size(U,1),1);
  else
      alpha = U(:, 2);
      U = U(:,1);
  end

  if length(U)==m.Ne
      U = U;
  elseif length(U)==m.Ne*m.Nq^2
    U = m.Qm\U;  % Least-Squares fit on the nodes for depiction
    alpha = alpha(1:m.Nn);
  end
  if length(alpha)==m.Ne
      alpha = alpha;
  elseif length(alpha)==m.Ne*m.Nq^2
    alpha = m.Qm\alpha;  % Least-Squares fit on the nodes for depiction
  end

  % Triangular Elements
  for e=1:m.Ne_Tri
    ei = m.Tri(e, 1);
    if any(e_on==ei)
      V = m.Nds(m.Tri(e, 2:end), :);
      if length(U)==m.Nn
	fill(V(:,1), V(:,2), U(m.Tri(e, 2:end)), 'EdgeColor', 'k', 'FaceAlpha', alpha(ei)); hold on
      elseif length(U)==m.Ne
	fill(V(:,1), V(:,2), U(ei), 'EdgeColor', 'k', 'FaceAlpha', alpha(ei)); hold on
      elseif length(U)==1
	fill(V(:,1), V(:,2), U, 'EdgeColor', 'k', 'FaceAlpha', alpha); hold on
      end
    end
  end

  % Quadrilateral Elements
  for e=1:m.Ne_Quad
    ei = m.Quad(e, 1);
    if any(e_on==ei)
      V = m.Nds(m.Quad(e, 2:end), :);
      if length(U)==m.Nn
	fill(V(:,1), V(:,2), U(m.Quad(e, 2:end)), 'EdgeColor', 'k', 'FaceAlpha', alpha(ei)); hold on
      elseif length(U)==m.Ne
	fill(V(:,1), V(:,2), U(ei), 'EdgeColor', 'k', 'FaceAlpha', alpha(ei)); hold on
      elseif length(U)==1
	fill(V(:,1), V(:,2), U, 'EdgeColor', 'k', 'FaceAlpha', alpha); hold on
      end
    end
  end
end
