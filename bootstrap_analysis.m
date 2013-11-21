function bootstrap_analysis(pr,n,a,b,s,q,t)
	%
	% Bootstrap Analysis
	%	Comparacion de metodos bootstrap. Se comparan intervalos
	%	de confianza con base en el metodo de bootstrap parametrico, 
	%	bootstrap no parametrico y cantidad pivotal. Para cada metodo
	%	se evaluan tanto la media como la media recortada para
	%	constrastar los resultados de manera mas adecuada.
	%	
	%	Parametros:
	%		pr	:= vector de parametros (default normal)
	%		n	:= tamano de la muestra inicial
	%		a	:= alpha := Intervalo de (1-a)100% de confianza
	%		b	:= numero de muestras simuladas en cada bootstrap
	%		s	:= numero de simulacinoes bootstrap efectuadas
	%		q	:= cantidad de analisis efectuados (muestras diferentes)
	%		t	:= parametro para la media recortada
	%
	close all;
	mu = pr(1);
	sigma = pr(2);
	mnum = 0;
	mm = true;
	for i = 1:4*q
		% Preparacion
		if mod(i-1,4) == 0 || mod(i-1,4) == 2
			% Necesario reiniciar vectores para los
			% nuevos calculos.
			sb = zeros(2,s);
			sb_2 = zeros(2,s);
			sbr = zeros(2,s);
			piv = zeros(1,2);
			pivr = zeros(1,2);
		end
		
		if mod(i-1,4) == 0
			% Empezando o ya pasaron los cuatro metodos:
			% necesario generar muestra nueva.
			mnum = mnum + 1;
			m = normrnd(mu,sigma,1,n);
			xbr = trimmean(m,t);
			xb = mean(m);
			for j = 1:s
				[sb(1,j),sb(2,j)] = bootstrap_p(pr,xb,n,a,b);
				[sbr(1,j),sbr(2,j)] = bootstrap_p(pr,xbr,n,a,b);
			end
		elseif mod(i-1,4) == 2
			% Ya pasaron los BPs, rehacer con BNP con la
			% misma muestra que los BPs.
			for j = 1:s
				[sb(1,j),sb(2,j)] = bootstrap_np(m,xb,n,a,b);
				[sbr(1,j),sbr(2,j)] = bootstrap_np(m,xbr,n,a,b);
			end
		end
		% IC Pivotal
		[piv(1),piv(2)] = pivotal(pr,xb,n,a);
		[pivr(1),pivr(2)] = pivotal(pr,xbr,n,a);
		% Generar sub-grafica
		subplot(2*q,2,i);
		hold on; grid on;
		plot([0,s],[mu,mu],'r');
		if mm == true
			% Grafica con media muestral
			plot([0,s],[xb,xb],'b--');
			p = barh((piv(2)+piv(1))/2,s,(piv(2)-piv(1)));
			e = errorbar(1:s, (sb(1,:)+sb(2,:))/2, ...
				(sb(1,:)+sb(2,:))/2 - sb(1,:), ...
				sb(2,:) - (sb(1,:)+sb(2,:))/2, 'o-');
			gi = min(min(min(sb(1,:)),piv(1)),mu)-.1;
			gs = max(max(max(sb(2,:)),piv(2)),mu)+.1;
			mm = false;
		else
			% Grafica con media muestral recortada
			plot([0,s],[xbr,xbr],'--');
			p = barh((pivr(2)+pivr(1))/2,s,(pivr(2)-pivr(1)));
			e = errorbar(1:s, (sbr(1,:)+sbr(2,:))/2, ...
				(sbr(1,:)+sbr(2,:))/2 - sbr(1,:), ...
				sbr(2,:) - (sbr(1,:)+sbr(2,:))/2, 'o-');
			gi = min(min(min(sbr(1,:)),pivr(1)),mu)-.1;
			gs = max(max(max(sbr(2,:)),pivr(2)),mu)+.1;
			mm = true;
		end
		if mod(i-1,4) == 0
			tl = title('B-P MM');
		elseif mod(i-1,4) == 1
			str = sprintf('B-P MM-R-%d',t);
			tl = title(str);
		elseif mod(i-1,4) == 2
			tl = title('B-NP MM');
		elseif mod(i-1,4) == 3
			str = sprintf('B-NP MM-R-%d',t);
			tl = title(str);
		end
		% Ajustes esteticos de graficas
		ch = get(p,'child');
		set(ch,'facea',.05);
		set(e, ...
			'LineStyle','none', ...
			'MarkerSize',3, ...
			'Color',[.3 .3 .3], ...
			'MarkerEdgeColor',[.2 .2 .2], ...
			'MarkerFaceColor',[.7 .7 .7]);
		errorbar_tick(e,150);
		set(tl,'FontSize',15);
		axis([-1, s+1, gi, gs]);
	end
end

function [ci,cs] = bootstrap_np(m,xb,n,a,b)
	m = sort(m);
	B = unidrnd(n,b,n);
	for i=1:b
		for j=1:n
			B(i,j) = m(B(i,j));
		end
	end
	y = zeros(b,1);
	for i=1:b
		y(i) = mean(B(i,:));
	end
	y = sort(y);
 	ci = 2*xb-y(max(1,round((1-a/2)*b)));
 	cs = 2*xb-y(max(1,round(a/2*b)));
end

function [ci,cs] = bootstrap_p(pr,xb,n,a,b)
	sigma = pr(2);
	B = normrnd(xb,sigma,b,n);
	y = zeros(b,1);
	for i=1:b
		y(i) = mean(B(i,:));
	end
	y = sort(y);
 	ci = 2*xb-y(max(1,round((1-a/2)*b)));
 	cs = 2*xb-y(max(1,round(a/2*b)));
end

function [ci,cs] = pivotal(pr,xb,n,a)
	sigma = pr(2);
	ci = xb+norminv(a/2,0,1)*sqrt(sigma)/sqrt(n);
	cs = xb-norminv(a/2,0,1)*sqrt(sigma)/sqrt(n);
end

function errorbar_tick(h,w,xtype)
	flagtype = get(h,'type');
	if nargin==1
		w = 80;
	end
	if nargin==2
	   xtype = 'ratio';
	end
	if ~strcmpi(xtype,'units')
		dx = diff(get(gca,'XLim'));
		w = dx/w;
	end
	if strcmpi(flagtype,'hggroup')
		hh=get(h,'children');
		x = get(hh(2),'xdata');
		x(4:9:end) = x(1:9:end)-w/2;
		x(7:9:end) = x(1:9:end)-w/2;
		x(5:9:end) = x(1:9:end)+w/2;
		x(8:9:end) = x(1:9:end)+w/2;
		set(hh(2),'xdata',x(:))
	else
		x = get(h(1),'xdata');
		x(4:9:end) = x(1:9:end)-w/2;
		x(7:9:end) = x(1:9:end)-w/2;
		x(5:9:end) = x(1:9:end)+w/2;
		x(8:9:end) = x(1:9:end)+w/2;
		set(h(1),'xdata',x(:))
	end
end
