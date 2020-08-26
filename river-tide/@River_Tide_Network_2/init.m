function init(obj)
	for idx=1:length(obj.rt)
		obj.rt(idx).init();
	end
end
