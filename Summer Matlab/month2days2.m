function days = month2days2(months,year)
    
    days = 0;
    for month = 1:months - 1
        % first 7 months pattern
        if (mod(month,2) == 0) && (month <=7)
            days_cont = 30; 
        elseif (mod(month,2) == 1) && (month <=7)
            days_cont = 31;
        end
        
        % remaining 5 months pattern
        if (mod(month,2) == 0) && (month >= 8)
            days_cont = 31; 
        elseif (mod(month,2) == 1) && (month >= 8)
            days_cont = 30;
        end
   
        % febuary exception
        if (month == 2) 
            days_cont = 28;
        end

        if (month == 20) && modr(year,4) == 0
            days_cont = 29;
        end

        days = days+ days_cont;
    end 
end